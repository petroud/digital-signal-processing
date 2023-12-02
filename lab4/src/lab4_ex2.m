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


%% Question (c)

Fs_C = 100; %Sampling Frequency
Ts_B = 1/Fs_B; 
