close all;
clear all;
clc;

%% Exercise 1

%Set constants of the low-pass filter
omega_c = 0.4*pi;  
Fc = omega_c/(2*pi); 
Fs = 100;
N = 21;

f_axis = @(w,Fs) 0:Fs/(2*length(w)):Fs/2-Fs/(2*length(w));


%Normalize the cutoff frequency
omega_n = Fc/(Fs/2);

%Create the filters. fir1() designs an N'th order lowpass FIR digital filter
% and returns the filter coefficients 
hFilter = fir1(N-1, omega_n, hamming(N));
rFilter = fir1(N-1, omega_n, rectwin(N));

%Get the frequency response at the N-point freq vector
[h1,w1] = freqz(hFilter,N);
[h2,w2] = freqz(rFilter,N);

%Fix the frequency axis 
hFreq = f_axis(w1,Fs);
rFreq = f_axis(w2,Fs);

%Plot the freq response
figure;
plot(hFreq, abs(h1));
hold on;
plot(rFreq, abs(h2));
xlabel('F(Hz)');
ylabel('Magnitude');
legend('Hamming Filter', 'Rectangular Filter');
title('Frequency Response of Hamming & Rectangular Filter');


