%% 1A

Fs = 8000;
T = 1/Fs;
StopTime = 1;
t = (0:T:StopTime-T)';
s = sin(2*pi*60*t);

what = sin(2*pi*2000*t);

a = 0;
r = .98;
u = .001;

% Transfer function coefficients
b0 = 1;
b1 = a;
b2 = 1;

a0 = 1;
a1 = a*r;
a2 = r^2;

b = [b0, b1, b2];
a = [a0, a1, a2];

[H, w] = freqz(b, a, 8000, 'whole', Fs);

% fvtool(x,1);

figure(1);
subplot(2, 1, 1);
plot(linspace(-pi, pi, length(s)), (fftshift(abs(fft(s+what)))));
title('Magnitude Response of s+w');
xlabel('rad/samp');
ylabel('Magnitude (dB)');

subplot(2, 1, 2);
plot(linspace(-pi, pi, length(s)), fftshift(angle(fft(s+what))));
title('Phase Response of s+w');
xlabel('rad/samp');
ylabel('Phase (Radians)');

figure(2);
plot(linspace(-pi, pi, length(H)), 20*log10(fftshift(abs(H))));
title('Magnitude Response of H');
xlabel('rad/samp');
ylabel('Magnitude (dB)');

shat=conv(s+what,ifft(H));

figure(3);
plot(linspace(-pi, pi, length(shat)), (fftshift(abs(fft(shat)))));
title('Magnitude Response of shat');
xlabel('rad/samp');
ylabel('Magnitude (dB)');
