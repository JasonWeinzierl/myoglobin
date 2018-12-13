clear
% change filename below as appropriate
rawdist=importdata('mb-ws.FEdist.dat');
figure, plot(rawdist), title('Fe displacement over time')
rawdistacorr=xcorr(rawdist);
poscorr=rawdistacorr(2001:4001);

%%%% uncomment the following to see frequency analysis of  position autocorrelation

% figure, plot(poscorr), title('position autocorrelation')
% y=poscorr;
% fbin=10000;
% Y = fft(y,fbin);
% Pyy = Y.* conj(Y) / fbin;
% f = 2*pi*(0:(fbin/2))/fbin;
% figure
% plot(f,Pyy(1:(fbin/2+1)))
% axis([0 0.120 0 5])
% title('Frequency content of position autocorrelation function')
% [C,I]=max(Pyy);, peakfreq=f(I)

units='angstroms,picoseconds'
gamma=3e-6
kray=7.3

timing=-2000:2000;
exphalfgammat=exp(-0.5*gamma*abs(timing));

distmax=max(rawdist)
distmin=min(rawdist)
scale=(min(rawdist):(max(rawdist)-min(rawdist))/99:max(rawdist));
histdist=hist(rawdist,100);
sig=std(rawdist)
a=mean(rawdist)
gauss=(9.3*exp(-1*(scale-a).^2/(2*sig^2))/sqrt(2*pi*sig^2));
figure, plot(scale,gauss/2001,scale,histdist/2001), title('Fe displacement distribution & gaussian fit')

kTinpNA=41.4
f=kTinpNA/sig^2
poscorr=rawdistacorr(2001:4001);
postime=timing(2001:4001);

% HARDCODED S,Q - use freq analysis as a starting point for S, fit by sight
% change as appropriate
%S=0.0044;, Q=0.365;
B = 0.0121;
S=0.0035;, Q=S/B;
for m=1:2001
  exp((m-1)*B/-2)*(cos((m-1)*Q*B)+(B/(2*Q*B))*sin((m-1)*Q*B));
  udamp(m)=ans;
  exp((m-1)*B/-2)*(cosh((m-1)*Q*B)+(B/(2*Q*B))*sinh((m-1)*Q*B));
  odamp(m)=ans;
  exp((m-1)*B/-2)*(1+B*(m-1)/2);
  nodamp(m)=ans;
end 
figure,plot(postime,poscorr,postime,odamp), axis([1 2000 -1 2]),title('position autocorrelation and fit')
w=sqrt(B^2/4 - S^2)
D=kTinpNA*w*w/(f*B)
alpha=S^2/B

fbin=40000;
for t=1:4001
        exp(-1*gamma*abs(t-2001)/2 - kray*(D/alpha)*(1-exp(-1*alpha*abs(t-2001))));
        tofft(t)=ans;
end
I = 0.015*fft(tofft,fbin)/2;
Pfft = I.* conj(I) / fbin;
mossfreq = 2*pi*(0:(fbin/2))/fbin;
negPfft = flipud(Pfft);
negmossfreq=-1*flipud(mossfreq);
figure
plot(mossfreq/gamma,Pfft(1:(fbin/2+1)),negmossfreq/gamma,negPfft(1:(fbin/2+1)))
title ('lineshape I(\omega), \omega in units of \Gamma')
axis ([-0.007/gamma 0.007/gamma -0.005 0.025])
