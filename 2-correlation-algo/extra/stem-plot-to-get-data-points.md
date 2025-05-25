
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

t_input = 0:T_input:0.1; % Time vector (0.1 seconds)

omega_input = 2 * pi * f0_input; % Angular frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Works with both data points and formula

orthogonal_function = sin(omega_input * t_input ); % Input waveform

t_input = (0:length(orthogonal_function)-1) * T_input; % Time vector

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots the ORIGINAL Signal and SAMPLES, BEFORE FILTERING

% Top subplot: Original continuous signal

figure;

subplot(1,1,1);

stem(t_input, orthogonal_function, 'r', 'filled');

title('Sampled Signal (Stems)');

xlabel('Time (s)');

ylabel('Sample Value');

grid on;
```