
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Parameters

fs_input = 720; % Sampling frequency (Hz)

T_input = 1 / fs_input; % Sampling period (s)

f0_input = 60; % Signal frequency (Hz)

samples = fs_input / f0_input;

half_samples = samples / 6;

t_input = T_input(1):T_input:half_samples / f0_input; % Time vector (duration = 1 cycle / 6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ORTHOGONAL SIGNAL. CHANGE THIS LINE TO OBTAIN DATA POINTS!

orthogonal_function = sin(2 * pi * 60 * t_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mapped symbolic signal

x_sign = sign(orthogonal_function); % -1, 0, or 1

% Even and Odd parts of mapped signal

x_sign_rev = fliplr(x_sign);

x_even = 0.5 * (x_sign + x_sign_rev);

x_odd = 0.5 * (x_sign - x_sign_rev);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Print x_sign values to console

fprintf('\nx_sign = [');

fprintf('%d, ', x_sign(1:samples-1));

fprintf('%d]\n', x_sign(samples));

% Plotting

figure;

stem(t_input, x_sign, 'k', 'filled');

title('Even Odd Functions');

xlabel('Time (s)');

ylabel('Sign Value');

ylim([-1.5 1.5]);

grid on;
```