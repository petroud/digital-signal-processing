clear all;
close all;
clc;

%% Exercise 1.A


%Set a space of 20 discrete point
n = -10:10;

x1 = ((n>=-4) - 2.*(n>=6));

figure(1)
subplot(4,1,1)
stem(n, x1);
title('\{x_1\}');
xlabel('n');
ylabel('x_1[n]');

x2 = 2*exp(n).*(n>=0);
subplot(4,1,2)
stem(n, x2);
title('\{x_2\}');
xlabel('n');
ylabel('x_2[n]');

%%Convolution

lenX1 = length(x1);
lenX2 = length(x2);

%Calculating the length of the convolved signal to be generated
nconv = 2*n(1): 2*n(end);
lenNconv = length(nconv);

%Add zero padding to the convolved signal
convSignal = zeros(1, lenNconv);

%%Implementing manual convolution using the sum rule as stated in the
%%exercise
for i=1:lenNconv
  for j = 0:lenX1-1
    if ( (i-j) > 0 && (i-j)<=lenX2 )
        convSignal(i) = convSignal(i) + x1(j+1).*x2(i-j);
    end
  end
end

subplot(4,1,3)
stem(nconv, convSignal)
title('Custom Convolution of \{x_1\}, \{x_2\}');
xlabel('n');
ylabel('Convolved Singal');

subplot(4,1,4)
stem(nconv, conv(x1, x2))
title('Convolution of \{x_1\}, \{x_2\} using [conv] function ');
xlabel('n');
ylabel('Convolved Singal');




%% Exercise 1.B

%Create 2 new discrete signals
x1_b = 4.*n;
x2_b = exp(n);


figure(2);
subplot(2,1,1)
stem(n,x1_b)
%x1=4n,n<=|10|
title('\{x_1\}*');
xlabel('n');
ylabel('x_1[n]*');

subplot(2,1,2)
stem(n,x2_b)
%x2=e^n n<=|10|
title('\{x_2\}*');
xlabel('n');
ylabel('x_2[n]*');


% Convolution using [conv]
z = conv(x1_b, x2_b);
nconv = 2 * n(1) : 2 * n(end);
figure(3)
subplot(2,1,1)
stem(nconv, z)
title('Convolution of \{x_1\}* and \{x_2\}*');
xlabel('n');
ylabel('z[n]');

% Multiplication in frequency field

% First, perform fourier transform to change signals into frequency field
X1_b = fft(x1_b, length(nconv));
X2_b = fft(x2_b, length(nconv));

% Multiply frequency fields
M = X1_b .* X2_b;

% Perform invert Fourier Transform to prove the lemma
z_f = ifft(M);

% Plotting the derived z[n]
subplot(2,1,2)
stem(nconv, z_f)
title('Reverse fourier transform of multiplied fourier of \{x_1\}* and \{x_2\}*');
ylabel('z[n]=F^{-1}\{Z[F]\}');
xlabel('n');