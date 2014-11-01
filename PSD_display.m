%written 12/2/13
%frequency 0<f<1/2; 1 corresponds to fs
%f0 is carrier frequency
%fb is the data frequency
%n is number of points in output array
% A is carrier amplitude

Eb=(10*log10((A^2/2)/(fb*fs)))*ones(size(y));
N0=(10*log10(2*sd^2/fs))*ones(size(y)); % SD of randn is 1
figure(1);
hold on;

axis([0 500 -70 0]);
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
title('Single-sided Periodogram and Welch PSD in dB/Hz versus frequency' );
hold off;