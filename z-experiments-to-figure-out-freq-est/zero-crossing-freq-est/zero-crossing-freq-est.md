
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 720; % Sampling frequency (Hz)

T_input = 1 / fs_input; % Sampling period (s)

f0_input = 60; % Signal frequency (Hz)

samples = fs_input/f0_input;

cycles = (samples/2)/f0_input;

t_input = 0:T_input:cycles; % Time vector (0.1 seconds)

Vm_input = 10; % Amplitude

omega_input = 2 * pi * f0_input; % Angular frequency

datapoints = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Works with both data points and formula

x = [714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699];

x = Vm_input * sin(omega_input * t_input + pi/18 ); % Input waveform

datapoints = false;

t_input = 0:T_input:((length(x)-1) * T_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots the ORIGINAL Signal and SAMPLES, BEFORE FILTERING

% Top subplot: Original continuous signal

figure;

subplot(2,1,1);

plot(t_input, x, 'b', 'LineWidth', 1);

title('Original Continuous Signal');

xlabel('Time (s)');

if(datapoints)

xlim([t_input(1), t_input(end)]);

xticks(0:T_input:t_input(end))

xtickformat('%.4f'); % shows more precise decimals

end

ylabel('Amplitude');

grid on;

% Bottom subplot: Sampled signal (stem)

subplot(2,1,2);

stem(t_input, x, 'r', 'filled');

title('Sampled Signal (Stems)');

xlabel('Time (s)');

if(datapoints)

xlim([t_input(1), t_input(end)]);

xticks(0:T_input:t_input(end))

xtickformat('%.4f'); % shows more precise decimals

end

ylabel('Sample Value');

grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

t_vector = (0:length(x)-1)/fs_input; % Time vector

% Sliding window parameters

window_duration = 0.05; % 50 ms

step_duration = 0.005; % 5 ms

window_size = round(window_duration * fs_input);

step_size = round(step_duration * fs_input);

num_steps = floor((length(x) - window_size) / step_size);

estimates = zeros(1, num_steps);

time_axis = zeros(1, num_steps);

fprintf("Zero Crossings:\n");

for k = 1:num_steps

start_idx = (k - 1) * step_size + 1;

end_idx = start_idx + window_size - 1;

x_win = x(start_idx:end_idx);

t_win = t_vector(start_idx:end_idx);

zero_crossings = [];

% zero-crossing detection and interpolation

for n = 2:length(x_win)

if x_win(n-1) * x_win(n) < 0

t_zero = t_win(n-1) - x_win(n-1)*(t_win(n) - t_win(n-1))/(x_win(n) - x_win(n-1));

zero_crossings(end+1) = t_zero;

% Only print for first few windows

if k <= 2

fprintf(" Window %d: t(%d) = %.6f s\n", k, length(zero_crossings), t_zero);

end

end

end

M = length(zero_crossings);

if M < 2

f_est = NaN; % Not enough crossings

else

t1 = zero_crossings(1);

tM = zero_crossings(end);

f_est = (M - 1) / (2 * (tM - t1));

end

estimates(k) = f_est;

time_axis(k) = t_win(floor(end/2)); % Center of window

end

fprintf("\nEstimated Frequency Range: %.2f Hz to %.2f Hz\n", min(estimates,[],'omitnan'), max(estimates,[],'omitnan'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plotting

figure;

plot(time_axis, estimates, 'b');

xlabel('Time (s)');

ylabel('Estimate (Hz)');

title('Zero Crossing Frequency Estimation');

grid on;
```




### Decaying DC as Input Signal:

![](20250601003141.png)
![](20250601003055.png)
![](20250601003037.png)
![](20250601003121.png)