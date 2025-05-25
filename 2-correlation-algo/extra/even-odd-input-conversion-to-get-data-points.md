
```
clc;

clear;

close all;

% Parameters

fs = 720;

T = 1 / fs;

t = T(1):T:0.02;

% Original signal (ONLY change this line)

x = sin(2 * pi * 60 * t);

% Mapped symbolic signal

x_sign = sign(x); % Will be -1, 0, or 1 based on the value of x

% Even and Odd parts of mapped signal

x_sign_rev = fliplr(x_sign);

x_even = 0.5 * (x_sign + x_sign_rev);

x_odd = 0.5 * (x_sign - x_sign_rev);

% Plotting

figure;

subplot(1,1,1);

stem(t, x_sign, 'k', 'filled');

title('Symbolic Signal (sign of x(t))');

xlabel('Time (s)');

ylabel('Sign Value');

ylim([-1.5 1.5]);

grid on;
```