close all
clear all
clc

%% Exercise 2

%Set time vector 0 to 500ms
t = 0:0.001:0.5;

%Set the required signal
y = 5.* cos(24.* pi .*t) - 2 .* sin(1.5.* pi .*t);


i=1;
for Ts=[1/48 1/24 1/12]
    t_sample=0:Ts:0.5;
    y_s = 5.* cos(24.* pi .*t_sample) - 2 .* sin(1.5.* pi .*t_sample);
    
    subplot(3,1,i);
    plot(t, y);
    grid on;
    hold on;
    stem(t_sample, y_s);
    
    hold off;
    title(['Sampling with T_s= ',num2str(Ts),' s']);
    xlabel('Time in seconds');
    ylabel('Amplitude');

    i=i+1;
end


A = 52;
Ts = 1/A;

t_sample = 0:Ts:0.5;

y_s = 5.* cos(24.* pi .*t_sample) - 2 .* sin(1.5.* pi .*t_sample);

figure(2)
plot(t, y)
title(['Sampling with T_s= 1/52 s']);
xlabel('Time in seconds');
ylabel('Amplitude');
grid on;
hold on;
stem(t_sample, y_s)
legend('y(t)', 'y_n[n] (sampling of y(t))')
hold off;