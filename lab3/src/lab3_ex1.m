clear all;
close all;
clc;

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

    %Create the digital filter
    digital_filter = freqz(numZ, denomZ, sample_space, Fs);

    
    %Plot the two filters in logarithmic scale
    figure;
    plot(sample_space, 20*log10(abs(analog_filter)),'r--');
    hold on;
    plot(sample_space, 20*log10(abs(digital_filter)),'b');
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
    title(sprintf('Butterworth Low-Pass Filter for Attenuation=%d dB', Rs));
    legend('Analog Filter', 'Digital Filter');
    grid on;
    hold off;
end

