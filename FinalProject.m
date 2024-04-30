%%
r = .95;
u = .000001;

Fs = 200;
T = 1/Fs;
StopTime = 500;
t = (0:T:StopTime-T)';
s = sin(2*pi*50*t);

i = 50*sin(2*pi*11*t);
x = s+i;

figure(1);
plot(linspace(-pi, pi, length(x)), 20*log10(fftshift(abs(fft(x)))));
title('X');
xlabel('radians/sec');
ylabel('in dB');

x = s+i;
y = zeros(1,length(x));
e = zeros(1,length(x));
a = zeros(1,length(x));

for index = 3:length(x)

    e(index) = x(index) + a(index) .* x(index-1) + x(index-2);
    y(index) = e(index) - r.*a(index).*y(index-1) - (r^2).*y(index-2);


    a(index+1) = a(index) - u.*y(index).*x(index-1);
    if (a(index+1) > 2) || (a(index+1) < -2)
        a(index+1) = a(index);
    end
end

figure(2);
subplot(2, 1, 1);
plot(linspace(0, length(y), length(y)), y);
title('y');
xlabel('sample');

len = length(y);
perc = .2;
subplot(2, 1, 2);
plot(linspace(-pi, pi, length(y(len*(1-perc):len))), 20*log10(fftshift(abs(fft(y(len*(1-perc):len))))));
title('mag of y');
xlabel('rad/samp');
ylabel('dB');

figure(3);
plot(linspace(0, length(a), length(a)), a);
title('a');
xlabel('time');

b0 = 1;
b1 = a(end);
b2 = 1;
a0 = 1;
a1 = a(end)*r;
a2 = r^2;
b = [b0, b1, b2];
a3 = [a0, a1, a2];
[H, w] = freqz(b, a, 10000, 'whole', Fs);
figure(4);
plot(linspace(-pi, pi, length(H)), 20*log10(fftshift(abs(H))));
title('Magnitude Response of H');
xlabel('rad/samp');
ylabel('Magnitude (dB)');

%%

chp = 50*chirp(t,0,500,100);
s = sin(2*pi*9*t);

x = s+chp;
% pspectrum(chp,t,'spectrogram','TimeResolution',0.1, ...
%     'OverlapPercent',99,'Leakage',0.85)
u = .000001;
R=.98;

figure(5);
plot(linspace(-pi, pi, length(s)), 20*log10(fftshift(abs(fft(x)))));
title('x');
xlabel('sample');
ylabel('dB');

y = zeros(1,length(x));
e = zeros(1,length(x));
a = zeros(1,length(x));
a(1) = -1.99999;
a(2) = -1.9999;
a(3) = -1.9999;

for index = 3:length(x)
    e(index) = x(index) + a(index) .* x(index-1) + x(index-2);
    y(index) = e(index) - r.*a(index).*y(index-1) - (r^2).*y(index-2);

    a(index+1) = a(index) - u.*y(index).*x(index-1);
    if (a(index+1) > 2) || (a(index+1) < -2)
        a(index+1) = a(index);
    end
end

figure(6);
plot(linspace(-pi, pi, length(y)), 20*log10(fftshift(abs(fft(y)))));
title('mag of y2');
xlabel('rad/samp');
ylabel('dB');

figure(7);
plot(linspace(0, length(a), length(a)), a);
title('a');
xlabel('time');

%%

r = .94;
u = .000005;

Fs = 200;
T = 1/Fs;
StopTime = 500;
t = (0:T:StopTime-T)';
s = sin(2*pi*8*t);

i1 = 50*sin(2*pi*57*t);
i2 = 50*sin(2*pi*67*t);

x = s+i1+i2;

figure(1);
plot(linspace(-pi, pi, length(x)), 20*log10(fftshift(abs(fft(x)))));
title('x');
xlabel('sample');
ylabel('dB');

y = zeros(1,length(x));
temp = zeros(1,length(x));
a1 = zeros(1,length(x));
a2 = zeros(1,length(x));

for index = 3:length(x)
    e = x(index) + a1(index) .* x(index-1) + x(index-2);
    temp(index) = e - r.*a1(index).*temp(index-1) - (r^2).*temp(index-2);
    a1(index+1) = a1(index) - u.*temp(index).*x(index-1);

  
    if (a1(index+1) > 2) || (a1(index+1) < -2)
        a1(index+1) = a1(index);
    end
end

figure(7);
plot(linspace(-pi, pi, length(temp(length(temp)*(1-perc):length(temp)))), 20*log10(fftshift(abs(fft(temp(length(temp)*(1-perc):length(temp)))))));
title('mag of y');
xlabel('rad/samp');
ylabel('dB');

for index = 3:length(x)

    e = temp(index) + a2(index) .* temp(index-1) + temp(index-2);
    y(index) = e - r.*a2(index).*y(index-1) - (r^2).*y(index-2);
    
    a2(index+1) = a2(index) - u.*y(index).*x(index-1);
    if (a2(index+1) > 2) || (a2(index+1) < -2)
        a2(index+1) = a2(index);
    end
end
figure(8);
plot(linspace(-pi, pi, length(y(length(y)*(1-perc):length(y)))), 20*log10(fftshift(abs(fft(y(length(y)*(1-perc):length(y)))))));
title('mag of y');
xlabel('rad/samp');
ylabel('dB');

figure(3);
subplot(2, 1, 1);
plot(linspace(0, length(a1), length(a1)), a1);
title('a1');
xlabel('time');

subplot(2, 1, 2);
plot(linspace(0, length(a2), length(a2)), a2);
title('a2');
xlabel('time');
