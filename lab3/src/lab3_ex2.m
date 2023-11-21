clear all;
close all;
clc;

%Set the constants for the Chebyshev filtering process
Ts=0.2;   %Sampling period
Fs=1/Ts;   %Sampling frequency
omega_c=2; %Cuttoff frequency
Rp=3;      %Passband ripple
N=256;     %Number of samples


%Create the required derived constants
omega_c_norm = omega_c/(Fs*pi);  %Normalized cutoff frequency

%Create highpass Chebyshev filter with grade 2
[num_2, denom_2] = cheby1(2, Rp, omega_c_norm, 'high');
[H2, omega2] = freqz(num_2, denom_2, N);


%Create highpass Chebyshev filter with grade 16
[num_16, denom_16] = cheby1(16, Rp, omega_c_norm, 'high');
[H16, omega16] = freqz(num_16, denom_16, N);


f = linspace(0, Fs/2, N);
f = 2*f/Fs; %Normalize sample freq

%Plot the filters in logarithmic dB
plot(f, 20*log10(abs(H2)), 'r--');
hold on;
plot(f, 20*log10(abs(H16)), 'b');
grid on;
legend('2nd order','16th order');
xlabel('Frequency (rad/sample)');
ylabel('Magnitude (dB)');
title('Chebyshev highpass filter');
xlim([-0.5 1.5]);
ylim([-250 20]);




