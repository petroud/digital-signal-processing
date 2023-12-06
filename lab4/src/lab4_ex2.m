clear all;
close all;
clc;

% Exercise 2 


%% Creation of common constants for all three questions
omega_c = 0.5*pi;
Fc = omega_c/(2*pi);
Nsmall = 21; 
Nbig = 41;

Fs = 100;

omega_n = Fc/(Fs/2);

%Set the axis fixing speedfunction
f_axis = @(w,Fs) 0:Fs/(2*length(w)):Fs/2-Fs/(2*length(w));


%Create the required filters, common for all quetsions

%Hammings
hamm21 = fir1(Nsmall-1, omega_n, hamming(Nsmall));
hamm41 = fir1(Nbig-1, omega_n, hann(Nbig));

%Hammings: Get the frequency response
[h_hamm21, w_hamm21] = freqz(hamm21,Nsmall);
[h_hamm41, w_hamm41] = freqz(hamm41,Nbig);


%-------------------------------

%Hannings
hann21 = fir1(Nsmall-1, omega_n, hann(Nsmall));
hann41 = fir1(Nbig-1, omega_n, hann(Nbig));

%Hannings: Get the frequency response
[h_hann21, w_hann21] = freqz(hann21,Nsmall);
[h_hann41, w_hann41] = freqz(hann41,Nbig);



%% Question (a)

Fs_A = 100; %Sampling Frequency

%Hammings: Fix the frequency axis
hamm21_f = f_axis(w_hamm21, Fs_A);
hamm41_f = f_axis(w_hamm41, Fs_A);

%Hammings: Fix the frequency axis
hann21_f = f_axis(w_hann21, Fs_A);
hann41_f = f_axis(w_hann41, Fs_A);

%Plot the hammings
figure();
subplot(1,2,1);
plot(hamm21_f, abs(h_hamm21));
title('Hamming filter with N=21');
ylabel('Magnitude');
xlabel('F(Hz');
subplot(1,2,2);
plot(hamm41_f, abs(h_hamm41));
title('Hamming filter with N=41');
ylabel('Magnitude');
xlabel('F(Hz');


%Plot the hannings
figure();
subplot(1,2,1);
plot(hann21_f, abs(h_hann21));
title('Hanning filter with N=21');
ylabel('Magnitude');
xlabel('F(Hz');
subplot(1,2,2);
plot(hann41_f, abs(h_hann41));
title('Hanning filter with N=41');
ylabel('Magnitude');
xlabel('F(Hz');



%% Question (b)

Fs_B = 100; %Sampling Frequency
Ts_B = 1/Fs_B;

%Set the sample space
N  = 2048;
n = 0:N-1;

%Set the signal
x = sin(15*n*Ts_B) + 0.25*sin(200*n*Ts_B);


%Set the new frequency axis needed after the Fourier Transform and the
%applied sampling
f_axis = -Fs/2:Fs/N:Fs/2-Fs/N;

%Calculate the FT of x(t)
X = fftshift(fft(x));

%Filter the signal using the 4 filters
x_hamm21 = filter(hamm21, 1, x);
x_hamm41 = filter(hamm41, 1, x);
x_hann21 = filter(hann21, 1, x);
x_hann41 = filter(hann41, 1, x);

%Calculate FT of each filtered result to get the spectrum
X_mm21 = fftshift(fft(x_hamm21));
X_mm41 = fftshift(fft(x_hamm41));
X_nn21 = fftshift(fft(x_hann21));
X_nn41 = fftshift(fft(x_hann41));

figure()
subplot(3,1,1);
plot(f_axis, abs(X));
xlabel('F(Hz)');
ylabel('|X(F)|')
title('Spectrum of x(t)(X(F))');

subplot(3,1,2);
plot(f_axis, abs(X_mm21));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hamming with N=21');

subplot(3,1,3);
plot(f_axis, abs(X_mm41));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hamming with N=41');


figure()
subplot(3,1,1);
plot(f_axis, abs(X));
xlabel('F(Hz)');
ylabel('|X(F)|')
title('Spectrum of x(t)(X(F))');

subplot(3,1,2);
plot(f_axis, abs(X_nn21));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hanning with N=21');

subplot(3,1,3);
plot(f_axis, abs(X_nn41));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hanning with N=41');




%% Question (c)

Fs_C = 50; %Sampling Frequency
Ts_C = 1/Fs_C;

%Set the sample space
N  = 2048;
n = 0:N-1;

%Set the signal
x = sin(15*n*Ts_C) + 0.25*sin(200*n*Ts_C);


%Set the new frequency axis needed after the Fourier Transform and the
%applied sampling
f_axis = -Fs_C/2:Fs_C/N:Fs_C/2-Fs_C/N;

%Calculate the FT of x(t)
X = fftshift(fft(x));

%Filter the signal using the 4 filters
x_hamm21 = filter(hamm21, 1, x);
x_hamm41 = filter(hamm41, 1, x);
x_hann21 = filter(hann21, 1, x);
x_hann41 = filter(hann41, 1, x);

%Calculate FT of each filtered result to get the spectrum
X_mm21 = fftshift(fft(x_hamm21));
X_mm41 = fftshift(fft(x_hamm41));
X_nn21 = fftshift(fft(x_hann21));
X_nn41 = fftshift(fft(x_hann41));

figure()
subplot(3,1,1);
plot(f_axis, abs(X));
xlabel('F(Hz)');
ylabel('|X(F)|')
title('Spectrum of x(t)(X(F)) at F_s=50Hz');

subplot(3,1,2);
plot(f_axis, abs(X_mm21));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hamming with N=21');

subplot(3,1,3);
plot(f_axis, abs(X_mm41));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hamming with N=41');


figure()
subplot(3,1,1);
plot(f_axis, abs(X));
xlabel('F(Hz)');
ylabel('|X(F)|')
title('Spectrum of x(t)(X(F)) at F_s=50Hz');

subplot(3,1,2);
plot(f_axis, abs(X_nn21));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hanning with N=21');

subplot(3,1,3);
plot(f_axis, abs(X_nn41));
xlabel('F(Hz)');
ylabel('|X(F)|');
title('X(F) filtered using Hanning with N=41');

