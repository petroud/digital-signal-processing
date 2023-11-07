close all;
clear all;
clc;

%% Exercise 1 

% Question b

% Set the frequency 
fs=1;
Ts=1/fs;

% Model the Transfer Function H(z) using tf()
H = tf( [0, 0.2 ,0], [1, -0.7, -0.18], Ts)

% Find the zeroes and the poles of the TF by solving the eqs
zeroes = roots([0, 0.2, 0])
poles = roots([1, -0.7, -0.18])

% Set the zeroes and poles on the Img plane
figure(1);
zplane(zeroes, poles);
title('Zeroes and Poles of the H(z)')



% Question d

%Set the frequency interval vector
fval = -pi:pi/128:pi;

figure(2);
% Normal drawing of frequency response
freqz( [0, 0.2, 0], [1, -0.7, -0.18], fval)
title('Frequency response for [-\pi,\pi] with step = \pi/128')

figure(3);
% Drawing of frequency response without providing interval
freqz( [0, 0.2, 0], [1, -0.7, -0.18])
title('Frequency response for [-\pi,\pi] with step = \pi/128 without setting the interval')




%Question e

clear all;

% Set the frequency 
fs=1;
Ts=1/fs;

% Model the Transfer Function H(z) using tf()
H = tf( [0, 0.2, 0], [1, -1.7, 0.52, 0.18], Ts)

% Find the zeroes and the poles of the TF by solving the eqs
zeroes = roots([0, 0.2, 0])
poles = roots([1, -1.7, 0.52, 0.18])

% Set the zeroes and poles on the Img plane
figure(4);
zplane(zeroes, poles);
title('Zeroes and Poles of the H*(z)');

%Set the frequency interval vector
fval = -pi:pi/128:pi;

figure(5);
% Normal drawing of frequency response
freqz( [0, 0.2, 0], [1, -1.7, 0.52, 0.18], fval)
title('Frequency response for [-\pi,\pi] with step = \pi/128')





%% Exercise 2

clear all;
clc;

% Set the numerator and denominator details
num = [4.0 -3.5 0];
denom = [1.0 -2.5 1.0];


% Question A

% Find zeroes and poles of the TF 
[r, p] = residuez(num, denom)


syms z

%Model the transfer function using simple fractions
H = r(1) / (1-p(1)*z^-1)  +  r(2) / (1-p(2)*z^-1);

%Print pretty the tf in the command line
pretty(H)


% Question B

% Calculate the inverse Z-Transform of H(z)
invZ = iztrans(H)






