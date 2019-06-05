% load('easy.mat');
function [decoded_text, y] = Rx(r, rolloff, desired_user)
%MORIMOTORX Summary of this function goes here
%OUTPUTS
%decoded_text: decoded message of the desired user
%y: output of equalizer
%INPUTSf
%r: received signal vector
%rolloff: beta
%desired_user: 1-3
%Parameters
% r=[r; zeros(10000,1)]; %For indexing
fc = 300e3;
Fs = 850e3;
Ts = 1/Fs;
IF = 2e6; %Assigned intermediate frequency
T = 6.4e-6; %Nominal symbol period
beta = rolloff; %SRRC pulse shape rolloff factor (.1 to .3)
slotAlot = 204e3; %FDM user slot allotment
psWid = 8; %8 clock periods
Pre = 'A0Oh well whatever Nevermind';
datLen = 875; %Symbols per user
numUsers = 3;
frameLen = 2870;

pSize = .5*8; %Half the length of the srrc pulse (in symbols)
P = T/Ts; %Oversampling factor (samples per symbol)
t_off = 0;
%IF Carrier phase = filtered white noise
%IF offset = fixed, < .01%
%Timing offset. Fixed.
%Symbol period offset. Fixed. 
IntersymbolInt = 7; %Max delay spread = 7 symbols
sampPerFrame = 4760;

%Carrier Recovery and Downconversion
%Squaring Data
sqrtDat = r.^2;


%BPF
wpass = [230e3 270e3]/(Fs/2);
rp= bandpass(sqrtDat, wpass);
% fcenter=250e3;
% fl = 200;
% fspan = 200;
% ff=[0 (fcenter-500) (fcenter-10) (fcenter+10) (fcenter+500) (Fs/2)]/(Fs/2);
% fa = [0 0 1 1 0 0];
% h = firpm(fl, ff, fa);
% rp = filter(h, 1, sqrtDat);
% figure(100)
% plotspec(rp, Ts);

%Estimate Phase offset of BPF
fftrBPF=fft(rp);                     % spectrum of rBPF
[m,imax]=max(abs(fftrBPF(1:end/2))); % find freq of max peak
ssf=(0:length(rp))/(Ts*length(rp));  % frequency vector
freqS=ssf(imax);                      % freq at the peak
phasep=angle(fftrBPF(imax));         % phase at the peak
[IR,f]=freqz(rp,1,length(rp),1/Ts);   % freq response of filter
[mi,im]=min(abs(f-freqS));           % freq where peak occurs
phaseBPF=angle(IR(im));            % < of BPF at peak freq
phaseS=mod(phasep-phaseBPF,pi);       % estimated angle


%Dual PLL 
time=5; t=0:Ts:length(r)*Ts-Ts; % time vector
mu1=.01; mu2=.003;                  % algorithm stepsizes
f0=300e3;                            % assumed freq at receiver
lent=length(r); th1=zeros(1,lent);  % initialize estimates
th2=zeros(1,lent); carest=zeros(1,lent);
for k=1:lent-1                      % combine top PLL th1
  th1(k+1)=th1(k)-mu1*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+phaseS);   %Add phaseS?        
  th2(k+1)=th2(k)-mu2*rp(k)*sin(4*pi*f0*t(k)+2*th1(k)+2*th2(k)+phaseS);  %Add phaseS?
  % with bottom PLL th2 to form estimate of preprocessed signal
  carest(k)=cos(4*pi*f0*t(k)+2*th1(k)+2*th2(k));
end
th = th1+th2;


% figure(1);
% subplot(3,1,1), plot(t,th1)              % plot first theta
% title('output of first PLL')
% ylabel('\theta_1')
% subplot(3,1,2), plot(t,th2)              % plot second theta
% title('output of second PLL')
% ylabel('\theta_2')
% subplot(3,1,3), plot(rp-carest) % plot difference between estimate
%                                 % and preprocesssed carrier rp
% title('error between preprocessed carrier and estimated carrier')
% xlabel('time')
% ylabel('f_0 - f')

%Downconversion
modi = 2*cos(2*pi*fc*t + th);
newDat = r.*modi';
% figure(2);
% plotspec(newDat, Ts);

% %Low pass filter
% freqLPF = 1.25e5;
% n = 200;
% f11 = [0 (freqLPF - 20000) (freqLPF+20000) (Fs/2)]/(Fs/2);
% a = [1 1 0 0];
% b = firpm(n, f11, a);
% downrLPF = filter(b, 1, newDat);
% figure(3);
% plotspec(downrLPF, Ts);


%Matched Filter (SRRC)
s1 = srrc(pSize, beta, P, t_off);
% s1 = fliplr(s1);
convOut = conv(s1, newDat);

%Interpolation downsampler
%Output is in SYMBOLS
% run clock recovery algorithm
n=ceil(length(newDat)/P);  % number of data points
tnow=pSize*P+1; tau=0; xs=zeros(1,n);   % initialize variables
tausave=zeros(1,n); tausave(1)=tau; i=0;
mu=999000;                            % algorithm stepsize
% delta=5e-8;                          % time for derivative, 0.05
delta=1e-8;
while (tnow<length(convOut)-(pSize*P*2))&&(i<=n)            % run iteration
  i=i+1;
  xs(i)=interpsinc(convOut,tnow+tau,pSize, rolloff);   % interp at tnow+tau
  x_deltap=interpsinc(convOut,tnow+tau+delta,pSize);  % value to right
  x_deltam=interpsinc(convOut,tnow+tau-delta,pSize);  % value to left
  dx=x_deltap-x_deltam;             % numerical derivative
  tau=tau+mu*dx*xs(i);              % alg update (energy)
  tnow=tnow+P; tausave(i)=tau;      % save for plotting
end
figure(4);
subplot(2,1,1), plot(xs(1:i-2),'b.')    % plot constellation diagram
title('constellation diagram');
ylabel('estimated symbol values')
subplot(2,1,2), plot(tausave(1:i-2))               % plot trajectory of tau
ylabel('offset estimates'), xlabel('iterations')

%Correlation. Peaks are 245 symb. too late bc of correlation. 
pamF = letters2pam2(Pre);%NO LONGER CORRELATION?!?! so no need -245?

%CORRELATE
corrDat1 = conv(xs, fliplr(pamF));
% [pks,ind] = findpeaks((abs(corrDat1),'MinPeakHeight',1200);
[pks, ind] = findpeaks(abs(corrDat1), 'MinPeakHeight', 1200);

%Equalizer
j = 0; 
eqOut = [];
eqConv = [];
n=12; 
eqF= zeros(length(ind), n);
while j< length(ind)
    j = j+1;
    s = pamF;                   % ideal channel output
%     chOut = xs(ind(j)-245+1:ind(j)+244-245+1); %New channel each run. Corrupted Preamble. 
    chOut = xs(ind(j)+1-245:ind(j)+244+1-245); %New channel each run. Corrupted Preamble. 

    f=zeros(n,1);           % initialize equalizer at 0
    mu=.0027; delta=n/2;             % stepsize and delay delta
    err=[];
    for i=n+1:length(pamF)                 % iterate
      rr=chOut(i:-1:i-n+1)';         % vector of received signal
      e=pamF(i-delta)-rr'*f;        % calculate error
      f=f+mu*e*rr;               % update equalizer coefficients
      f;
    end
    eqF(j, :) = f.'; 
end

% % Output of equalizer
% j = 0;
% lastBlock = zeros(1, 2870);
% while j<length(ind)
%     j = j+1;
% %     v = [xs(ind(j)-245:ind(j)-245+2625+244)]; %need prev. frame for block conv.
%     v = [xs(ind(j)-245:ind(j)-245+2625+244)]; %need prev. frame for block conv.
% 
% %     eqConv = filter(eqF(j,:), 1, [lastBlock v]);
%     eqConv = conv(eqF(j,:), [lastBlock v]);
%     eqOut = [eqOut eqConv(2871:end)];
%     lastBlock = v;
% end
% Output of equalizer
j = 0;
lastBlock = zeros(1, 2870);
while j<length(ind)
    j = j+1;
%     v = [xs(ind(j)-245+1:ind(j)-245+2625+244+1)]; 
    v = [xs(ind(j)+1-245:ind(j)+2870-245)]; 
%     eqConv = filter(eqF(j,:), 1, [lastBlock v]);
    eqConv = conv(eqF(j,:), v);
    eqConv = eqConv(delta+1:2870+delta);
    eqOut = [eqOut eqConv];
end

%CORRELATE round 2 (rob)
corrDat2 = conv(eqOut, fliplr(pamF));
% [pks,ind] = findpeaks((abs(corrDat1),'MinPeakHeight',1200);
[pks, ind2] = findpeaks(abs(corrDat2), 'MinPeakHeight', 1200);


%Quantize
quantDat = quantalph(eqOut, [-3, -1, 1, 3]);
quantDat = quantDat.';
data = [];
k=0;
while k<length(ind2)
    k = k+1;
%     data = [data quantDat(ind(k)-245+1+245+(desired_user-1)*875:ind(k)-245+1+244+(desired_user)*875)];
    data = [data quantDat(ind2(k)+1+((desired_user-1)*875):ind2(k)+((desired_user)*875))]; %ADDED DELTA

end
message = pam2letters2(data)

decoded_text = message;
y = eqOut;
figure(5); plot(y, '.')
% % % 
end

