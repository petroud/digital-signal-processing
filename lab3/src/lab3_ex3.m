clear all;
close all;

%% Exercise 1

%Set Frequencies and Filter Details
Fs=10e3; %Sampling freq
Wp=2*pi*3e3; %Passband freq
Ws=2*pi*4e3; %Stopband freq

Rp=3; %Ripple


%Perform two series of calculations 1st for 30dB 2nd for 50dB
for Rs=[30,50]

    %Find butterwoth order
    [N,Wn] = buttord(Wp, Ws, Rp, Rs, 's');
    
    Rs, N

    %Find poles and zeroes of the analog filter
    [z, p, k] = buttap(N);

    %Construct TF using Zs and Ps;
    [num, denom] = zp2tf(z, p, k);


    %Transform the lowpass analog filter to lowpass analog filter
    %with cutoff freq = Wn
    [numS, denomS] = lp2lp(num, denom, Wn);

    %Set the sampling space by sampling frequency and number of samples
    sample_space = linspace(0, Fs/2, 2048);

    %Create the analog filter
    analog_filter = freqs(numS, denomS, 2*pi*sample_space);

    %Bilinear transformation with frequency prewarping
    %From s-domain TF to z-transform discrete equivalent
    [numZ, denomZ] = bilinear(numS, denomS, Fs);

    %Keeping those necessary for exercise 3
    if Rs == 30
        numZ_30 = numZ;
        denomZ_30 = denomZ;
    else
        numZ_50 = numZ;
        denomZ_50 = denomZ;
    end
    
end


%% Exercise 2

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


%% Question 3a
fs = 10^4;
Ts = 1/fs;
Nsamples = 500;
n = 0:Nsamples-1;
x = 1 + cos(1000 * n*Ts) + cos (16000 * n*Ts) + cos (30000 * n*Ts);

figure;
subplot(4,1,1)
plot(n,x)
title('500 Samples of x(t) at f_s=10kHz');

f_x = -fs/2:fs/Nsamples:fs/2 -fs/Nsamples; %Set axis for FT
F_x = fftshift(fft(x))*Ts; %Shift the lobes to the center and normalize them

subplot(4,1,2)
plot(f_x, abs(F_x))
ylabel('|X(F)|');
xlabel('F (Hz)');
title('Original Spectrum |X(F)| of x(t)');

%Filter the signal using the TFs of the butterworth digital filters
x_filtered30 = filter(numZ_30, denomZ_30, x);
x_filtered50 = filter(numZ_50, denomZ_50, x);

%Shift the lobes and normalize
x_filtered30_F = fftshift(fft(x_filtered30))*Ts;
x_filtered50_F = fftshift(fft(x_filtered50))*Ts;

subplot(4,1,3)
plot(f_x, abs(x_filtered30_F))
ylabel('|X(F)|');
xlabel('F (Hz)');
title('|X(F)| filtered using 30dB lowpass Butterworth digital filter');

subplot(4,1,4)
plot(f_x, abs(x_filtered50_F))
title('|X(F)| filtered using 50dB lowpass Butterworth digital filter');
ylabel('|X(F)|');
xlabel('F (Hz)');

%% Question 3b

Ts = 0.2;
fs = 1/Ts;
Nsamples = 500;
n = 0:Nsamples-1;

x = 1 + cos(1.5 * n*Ts) + cos(5 * n*Ts);

figure;
subplot(3, 1, 1)
plot(n, x)
title('500 Samples of x(t) at f_s=5Hz');

f_x = -fs/2:fs/Nsamples:fs/2 -fs/Nsamples; %Set axis for FT
x_spectrum = fftshift(fft(x))*Ts;

subplot(3,1,2)
plot(f_x, abs(x_spectrum));
ylabel('|X(F)|');
xlabel('F (Hz)');
title('Original Spectrum |X(F)| of x(t)');


%Filtering the signal using the 16th order Chebyshev highpass filter.
x_filtered3 = filter(num_16, denom_16, x); 
x_filtered3_f = fftshift(fft(x_filtered3))*Ts;
subplot(3,1,3)
plot(f_x, abs(x_filtered3_f));
title('|X(F)| filtered using 16th order highpass Chebyshev filter');
ylabel('|X(F)|');
xlabel('F (Hz)');
