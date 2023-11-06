%%%% Lab 1 %%%%

clear all;
close all;
%%Question 1A

n = -10:10;
x1 = ((n>=-2) - 2.*(n>=6));


figure(1)
stem(n, x1)
title('First signal of our choice');
xlabel('n');
ylabel('x1[n]');

x2 = n.*(n>=0);
figure(2)
stem(n, x2)
title('Second signal of our choice');
xlabel('n');
ylabel('x2[n]');

%%Convolution
lenX1 = length(x1);
lenX2 = length(x2);

%% nconv = n_x1(1) + n_x2(1) : n_x1(end) + n_x2)(end)
nconv = 2*n(1): 2*n(end);
lenNconv = length(nconv)

%%zero padding conv result
convSignal = zeros(1, lenNconv)

%%Simulating Sum as shown in lab pdf
for i=1:lenNconv
  for j = 0:lenX1-1
    if ( (i-j) > 0 && (i-j)<=lenX2 )
        convSignal(i) = convSignal(i) + x1(j+1).*x2(i-j);
    end
  end
end

figure(3)
stem(nconv, convSignal)
title('Custom Convolution of x1,x2');
xlabel('n');
ylabel('convSignal');


matlabConvSignal = conv(x1, x2);
figure(4)
stem(nconv, matlabConvSignal)
title('Convolution of x1,x2 using conv function ');
xlabel('n');
ylabel('matlabConvSignal[n]');

%%Question 1B
x1_b = 2.*n;
x2_b = 2 .*exp(n);


figure(5);
stem(n,x1_b)
title('x1_b[n]=2*n, n<=|10|');
xlabel('n');
ylabel('x1_b[n]');


figure(6);
stem(n,x2_b)
title('x2_b[n]=2*(e^n), n<=|10|');
xlabel('n');
ylabel('x2_b[n]');

%%%% Convolution
z = conv(x1_b, x2_b);
nconv = 2 * n(1) : 2 * n(end);
figure(7)
stem(nconv, z)
title('Convolution of x1_b and x2_b');
xlabel('n');
ylabel('z[n]');
%%%% Multiplication in frequency field

%%First, fourier transform to go into the frequency field
X1_b = fft(x1_b, length(nconv));
X2_b = fft(x2_b, length(nconv));

%%Then multiplication here
M = X1_b .* X2_b;

%%Finally, reverting to time field
z_f = ifft(M);

%%Plotting to check if it's the same signal( = same operation)
figure(8)
stem(nconv, z_f)
title('Reverse fourier transform of z[n]');
xlabel('z_f[n]');
ylabel('n');


%%Question 2

t = 0:0.001:0.5;
y = 5.* cos(24.* pi .*t) - 2 .* sin(1.5.* pi .*t);


%%%SubQuestion a
Ts = 1/48;

t_sampling = 0:Ts:0.5

y_n = 5.* cos(24.* pi .*t_sampling) - 2 .* sin(1.5.* pi .*t_sampling);

figure(9)
plot(t, y)
title('Same plot for y(t) and sampling(y_n[n]) for T_s=1/48');
xlabel('Time in seconds');
ylabel('Amplitude');
grid on;
hold on;
stem(t_sampling, y_n)
legend('y(t)', 'y_n[n] (sampling of y(t))')
hold off;


%%%SubQuestion b
Ts = 1/24;

t_sampling = 0:Ts:0.5

y_n = 5.* cos(24.* pi .*t_sampling) - 2 .* sin(1.5.* pi .*t_sampling);

figure(10)
plot(t, y)
title('Same plot for y(t) and sampling(y_n[n]) for T_s=1/24');
xlabel('Time in seconds');
ylabel('Amplitude');
grid on;
hold on;
stem(t_sampling, y_n)
legend('y(t)', 'y_n[n] (sampling of y(t))')
hold off;


%%%SubQuestion c
Ts = 1/12;

t_sampling = 0:Ts:0.5

y_n = 5.* cos(24.* pi .*t_sampling) - 2 .* sin(1.5.* pi .*t_sampling);

figure(11)
plot(t, y)
title('Same plot for y(t) and sampling(y_n[n]) for T_s=1/12');
xlabel('Time in seconds');
ylabel('Amplitude');
grid on;
hold on;
stem(t_sampling, y_n)
legend('y(t)', 'y_n[n] (sampling of y(t))')
hold off;



%%%SubQuestion d
A = 52;
Ts = 1/A;

t_sampling = 0:Ts:0.5

y_n = 5.* cos(24.* pi .*t_sampling) - 2 .* sin(1.5.* pi .*t_sampling);

figure(12)
plot(t, y)
title('Same plot for y(t) and sampling(y_n[n]) for T_s=1/52');
xlabel('Time in seconds');
ylabel('Amplitude');
grid on;
hold on;
stem(t_sampling, y_n)
legend('y(t)', 'y_n[n] (sampling of y(t))')
hold off;


%%%Question 3a

w=10.* cos(2.* pi.*20.* t) - 4.* sin(2.*pi.*40.* t + 5);

f1 = 20;
f2 = 40; %% fmax

fs = 500; %%fs >= 2fmax(= 80Hz)
Ts = 1/fs;
Nsamples = 128;

n = 0:Nsamples - 1;
w_sampling = 10.* cos(2*pi*f1*n*Ts) - 4.* sin(2*pi*f2*n*Ts + 5);


f_axis = -fs/2:fs/Nsamples:fs/2 - fs/Nsamples;
W = fftshift(fft(w_sampling));
figure(13)
stem(f_axis, W)
xlabel('Frequency(Hz)','Fontsize',12);
ylabel('W(F)','Fontsize',12);
title('FT of sampled W','FontSize',14);

%%%Question 3b
t = 0:0.01:0.5;
phi = 52;
fs = 8000;
Ts = 1/fs;
t_axis = 0:Ts:0.5;

N = 1024;
F = (-fs/2:fs/N:fs/2-fs/N);
count=1;
figure(14)

for f=(100:125:475)
  x = sin(2*pi*f*t_axis + phi);
  X=fftshift(fft(x,N)*Ts);
  subplot(4,1,count)
  count=count+1;
  plot(F,abs(X));
  title(['Spectrum with f = ' num2str(f)])
  xlabel('f(Hz)');
  ylabel('|X(F)|');
end;

count=1;
figure(15)
for f=(7525:125:7900)
  x = sin(2*pi*f *t_axis + phi);
  X=fftshift(fft(x,N)*Ts);
  subplot(4,1,count)
  count=count+1;
  plot(F,abs(X));
  title(['Spectrum with f = ' num2str(f)])
  xlabel('f(Hz)');
  ylabel('|X(F)|');
end;






