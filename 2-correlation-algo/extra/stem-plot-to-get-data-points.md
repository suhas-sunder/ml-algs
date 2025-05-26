
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 720; % Sampling frequency (Hz)

T_input = 1 / fs_input; % Sampling period (s)

f0_input = 60; % Signal frequency (Hz)

samples = fs_input / f0_input;

half_samples = samples / 6; % I made this more than half because we don't need to see more than 1-2 cycles for this plot to get our data samples

t_input = 0:T_input:half_samples / f0_input; % Time vector (duration = 1 cycle / 6)

omega_input = 2 * pi * f0_input; % Angular frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ORTHOGONAL SIGNAL. CHANGE THIS LINE TO OBTAIN DATA POINTS!

orthogonal_function = sin(omega_input * t_input); % Input waveform

t_input = (0:length(orthogonal_function)-1) * T_input; % Re-adjusted time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Display Sampled Values in Array Format

fprintf('\northogonal_function = [');

fprintf('%.4f, ', orthogonal_function(1:samples - 1));

fprintf('%.4f]\n', orthogonal_function(samples));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting the Signal as Stems

figure;

stem(t_input, orthogonal_function, 'r', 'filled');

title('Sampled Signal (Stems)');

xlabel('Time (s)');

ylabel('Sample Value');

grid on;
```