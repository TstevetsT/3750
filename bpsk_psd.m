function bpsk_psd
%written 12/2/13
%modified 4/8/93 to conform to matlab frequency definition,
%computes array of samples from a bpsk signal with amplitude 1
%frequency 0<f<1/2; 1 corresponds to fs
%f0 is carrier frequency
%fb is the data frequency
%n is number of points in output array

clear all
format compact;
fs=input('input sample frequency in Hz, R: ');
f0 =input(' Input digital carrier frequency, |f0| <= 1/2, R: ');
fb = input(' Input digital data frequency, |fb| <= |f0| <= 1/2, R: ');
n = input(' Input number of samples, n > 1, I: ');

A=sqrt(2);
sd=1;
t = 0:n-1;  %in samples

%	Form the carrier

x = A*cos(2*pi*f0*t);

%	Produce a random sequence of {0,1}
%rand('uniform');
j = 1;
for i=1:n;
 if fix((i-1)*fb) < fix(i*fb);
  r=rand ;
  if r <0.5;
     j=-1;
   else; 
     j=1;
  end;
 end;
 p(i)=j;
end;

randn;
y = p.*x;

figure;
Eb=(10*log10((A^2/2)/(fb*fs)))*ones(size(y));



hold on;
w=sd*randn(size(y));
z=y+w;
save ('bpskfil','f0','fb','n','y','p','z','w') ;


N0=(10*log10(2*sd^2/fs))*ones(size(y));

axis([0 500 -70 0]);
hold on;
ZZ=(1/fs)*fft(z,n);
ZD=10*log10(2*(fs/n)*(abs(ZZ).^2));

f=0:fs/n:(n/2-1)*fs/n;
plot(f,ZD(1:n/2),'y');
plot(Eb,'r');
plot(N0,'g');

grid on;

 h = spectrum.welch;    % Create a Welch spectral estimator. 
 Hpsd = psd(h,z,'fs',fs);             % Calculate the PSD 
             plot(Hpsd)  % Plot the PSD.
title('Single-sided Periodogram and Welch PSDin dB/Hz versus frequency' );
hold off;