clear all;
close all;
clc;


%% Exercise 3
t = 0:0.01:0.5;
x = 10.* cos(2.* pi.*20.* t) - 4.* sin(2.*pi.*40.* t + 5);


figure(1);
subplot(2,1,1)
plot(t,x);
xlabel('Time(s)');
ylabel('x(t)');
title('\{x(t)\}');


Fs = 200; % fs >= F_Nyquist = 2f_max = 80Hz
Ts = 1/Fs;

t_sample = 0:Ts:0.5; 
x_sample = 10.* cos(2*pi*20*t_sample) - 4.* sin(2*pi*40*t_sample + 5);

Nsamples = 128;
f_axis = [-Fs/2: Fs/Nsamples : Fs/2 - Fs/Nsamples];
X = fftshift(fft(x_sample, Nsamples)*Ts);

subplot(2,1,2);
plot(f_axis, X)
xlabel('Frequency(Hz)');
ylabel('X(F)');
title('FT of sampled X at 128 samples');


%% Question 3b

t = 0:0.01:0.5;
phi = 52;
Fs = 8000;
Ts = 1/Fs;
t_axis = 0:Ts:1;

N = 1024;
F = (-Fs/2:Fs/N:Fs/2-Fs/N);
count=1;
figure();
for f=(100:125:475)
  x = sin(2*pi*f*t_axis + phi);
  X=fftshift(fft(x,N)*Ts);
  subplot(4,1,count)
  plot(F,abs(X));
  title(['Spectrum for f = ' num2str(f) 'Hz'])
  xlabel('f(Hz)');
  ylabel('|X(F)|');
 
  count=count+1;
end;

%%
count=1;
figure()
for f=(7525:125:7900)
  x = sin(2*pi*f*t_axis + phi);
  X=fftshift(fft(x,N)*Ts);
  subplot(4,1,count)
  plot(F,abs(X));
  title(['Spectrum for f = ' num2str(f) 'Hz'])
  xlabel('f(Hz)');
  ylabel('|X(F)|');
 
  count=count+1;
end