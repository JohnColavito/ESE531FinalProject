%% 1A

N = 1012;
Fs = 8000;
T = 1/Fs;
StopTime = 1;
t = (0:T:StopTime-T)';
s = sin(2*pi*60*t);

what = sin(2*pi*2000*t);
length(what)

a = 1;
r = .98;
u = .05;

% Transfer function coefficients
b0 = 1;
b1 = a;
b2 = 1;

a0 = 1;
a1 = a*r;
a2 = r^2;

b = [b0, b1, b2];
a = [a0, a1, a2];

[H, w] = freqz(b, a, N, 'whole', Fs);

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
ylabel('Magnitude');

for c = 0:10
    
end

%%
r = .90;
u = .000001;

N = 1012;
T = 1/Fs;
StopTime = 10;
t = (0:T:StopTime-T)';
s = sin(2*pi*60*t);

i = 10*sin(2*pi*200*t);


figure(1);
subplot(2, 1, 1);
plot(linspace(-pi, pi, length(s)), s+i);
title('s');
xlabel('sample');

subplot(2, 1, 2);
plot(linspace(-pi, pi, length(s)), 20*log10(fftshift(abs(fft(s+i)))));
title('s+i');
xlabel('sample');

x = s+i;
y = zeros(1,length(x));
e = zeros(1,length(x));
a = zeros(1,length(x));

for index = 3:length(x)
    e(index) = x(index) + a(index) * x(index-1) + x(index-2);
    y(index) = e(index) - r*a(index)*y(index-1) - (r^2)*y(index-2);

    a(index+1) = a(index) - u*y(index)*x(index-1);
    if a(index+1) > 2 || a(index+1) < -2
        a(index+1) = 0;
    end
end

figure(2);
subplot(2, 1, 1);
plot(linspace(-pi, pi, length(y)), y);
title('Magnitude Response of s+i');
xlabel('rad/samp');
ylabel('Magnitude (dB)');

subplot(2, 1, 2);
plot(linspace(-pi, pi, length(y)), 20*log10(fftshift(abs(fft(y)))));
title('mag of y');
xlabel('rad/samp');
ylabel('dB');

figure(3);
plot(linspace(-pi, pi, length(a)), fftshift(abs(fft(a))));
title('a');
xlabel('time');

figure(4);
plot(linspace(0, length(a), length(a)), a);
title('a');
xlabel('time');