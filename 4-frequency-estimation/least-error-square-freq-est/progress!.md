
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 1930; % Sampling frequency (Hz)

T_input = 1 / fs_input; % Sampling period (s)

f0_input = 55; % Signal frequency (Hz)

t_input = 0:T_input:((fs_input/f0_input)/2)/f0_input; % Time vector (0.1 seconds)

Vm_input = 10; % Amplitude

omega_input = 2 * pi * f0_input; % Angular frequency

datapoints = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters FOR FILTER

%#ok<*UNRCH>

fs = fs_input; % Sampling frequency

f0 = 60; % Filter Frequency

T = 1 / fs; % Sampling period

omega = 2 * pi * f0; % Discrete angular frequency (radians/sample)

f_range = linspace(0, fs, 1000); % Frequency range to plot full range. This is only relevant later on when assigning z = exp(j omega T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data points

% x = [714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699];

x = Vm_input * sin(omega_input * t_input ); % Input waveform

datapoints = false;

estimated_freq = f0 * ones(1, length(t_input));

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

% LES filter generator

filter_choice = 0; % Default initialization

track_active_filter = []; % For error handeling for user selection of filter

% Decide which filter to ACTIVATE. Manually configured filters.

fundamental_filter = true;

second_harmonic_filter = true;

third_harmonic_filter = true;

fourth_harmonic_filter = true;

fifth_harmonic_filter = false;

sixth_harmonic_filter = false;

dc_filter = true;

% Add to tracker only if true

if fundamental_filter, track_active_filter(end + 1) = 1; end

if second_harmonic_filter, track_active_filter(end + 1) = 2; end

if third_harmonic_filter, track_active_filter(end + 1) = 3; end

if fourth_harmonic_filter, track_active_filter(end + 1) = 4; end

if fifth_harmonic_filter, track_active_filter(end + 1) = 5; end

if sixth_harmonic_filter, track_active_filter(end + 1) = 6; end

if dc_filter, track_active_filter(end + 1) = 6; end

%-----------------------------------------------------------

% Decide which filter to TARGET/APPLY (Only activate one)

filters = {"fundamental", "2ndHarmonic", "3rdHarmonic", "4thHarmonic", "5thHarmonic", "6thHarmonic", "DC"};

filter_choice = 1; % FUNDAMENAL

% filter_choice = 2; % 2nd Harmonic

% filter_choice = 3; % 3rd Harmonic

% filter_choice = 4; % 4th Harmonic

% filter_choice = 5; % 5th Harmonic

% filter_choice = 6; % 6th Harmonic

% filter_choice = 7; % DC

if ~ismember(filter_choice, track_active_filter)

error('YOU CHOSE A FILTER THAT YOU DID NOT ACTIVATE!');

end

target_filter = filters{filter_choice};

%-----------------------------------------------------------

sample_offset = 1; % Change this to whatever you need. Prof said it should be minimum 1.

% Initialize sample count

samples = 0 + sample_offset;

% Add 2 for each active filter.

if fundamental_filter, samples = samples + 2; end

if second_harmonic_filter, samples = samples + 2; end

if third_harmonic_filter, samples = samples + 2; end

if fourth_harmonic_filter, samples = samples + 2; end

if fifth_harmonic_filter, samples = samples + 2; end

if sixth_harmonic_filter, samples = samples + 2; end

if dc_filter, samples = samples + 2; end

window = samples; % In case you want custom window size, can change it here manually

nested_matrix_pinv_A = cell(1, length(t_input));

matrix_A = [];

real_values = []; % Real part of filter. H(z) Vp Cos(Theta)

imaginary_values = []; % Imaginary part of filter. H(z) Vp Sin(Theta)

for k = 1:length(t_input)

z_power = -1 * (window - 1)/2; % Assumes window will always be an odd number so that the center becomes V0, then even amouns before and after.

for n = 1:window

array_of_equations = []; % Reset at the beginning of each loop

if fundamental_filter

array_of_equations(end + 1) = sin(2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(2*pi*estimated_freq(k)* T * z_power);

end

if second_harmonic_filter

array_of_equations(end + 1) = sin(2 *2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(2 *2*pi*estimated_freq(k)* T * z_power);

end

if third_harmonic_filter

array_of_equations(end + 1) = sin(3 *2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(3 *2*pi*estimated_freq(k)* T * z_power);

end

if fourth_harmonic_filter

array_of_equations(end + 1) = sin(4 *2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(4 *2*pi*estimated_freq(k)* T * z_power);

end

if fifth_harmonic_filter

array_of_equations(end + 1) = sin(5 *2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(5 *2*pi*estimated_freq(k)* T * z_power);

end

if sixth_harmonic_filter

array_of_equations(end + 1) = sin(6 *2*pi*estimated_freq(k)* T * z_power);

array_of_equations(end + 1) = cos(6 *2*pi*estimated_freq(k)* T * z_power);

end

if dc_filter

array_of_equations(end + 1) = 1;

array_of_equations(end + 1) = 1 * z_power;

end

% Insert nth row into matrix_A at row `n`

matrix_A(n, :) = array_of_equations;

z_power = z_power + 1;

end

target_array_location = 1;

% Perform left pseudo inverse on matrix A

A_left_pinv = pinv(matrix_A);

nested_matrix_pinv_A{k} = A_left_pinv;

end

%-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

window_size = samples; % Window size is same as sample size

V_real = zeros(1, window_size);

V_imaginary = zeros(1, window_size);

% Allocate arrays to store angles and magnitude values

phase_angle_deg = zeros(1, length(t_input));

mag = zeros(1, length(t_input));

% Apply the 3-sample phasor magnitude and angle estimator

% This is where we actually take 3 SAMPLES and APPLY THE FILTER

x_buffer = zeros(1, window_size); % sliding buffer for x

for n = 1:length(t_input)

target_pinv_A = nested_matrix_pinv_A{n} ./ 100000; % Use curly braces to get the matrix out

real_values = target_pinv_A(target_array_location, :);

imaginary_values = target_pinv_A(target_array_location + 1, :);

disp(target_pinv_A)

disp("real")

disp(real_values )

disp("im")

disp(imaginary_values )

% Slide buffer: drop oldest, append new sample

x_buffer = [x_buffer(2:end), x(n)]; % This takes all elements of x_buffer except the first one, destructures the array, then we add the result of x(n) at the end

% Apply weights to current buffer

V_real = real_values .* x_buffer;

V_imaginary = imaginary_values .* x_buffer;

imaginary_part_Vp_cos_theta = sum(V_imaginary);

real_part_Vp_sin_theta = sum(V_real);

phase_angle_deg(n) = atan(imaginary_part_Vp_cos_theta/ real_part_Vp_sin_theta) * 180 / pi;

mag(n) = sqrt(imaginary_part_Vp_cos_theta^2 + real_part_Vp_sin_theta^2);

if(n == window_size)

% Insert real and imaginary filter values into columns of table

array_for_table(1:window_size, 2) = round(real_values(1:window_size), 4);

array_for_table(1:window_size, 5) = round(imaginary_values(1:window_size), 4);

% Insert calculated values into columns for table after calculation is done.

array_for_table(1:window_size, 3) = round(V_real(1:window_size));

array_for_table(1:window_size, 6) = round(V_imaginary(1:window_size));

array_for_table(1:window_size, 4) = round(cumsum(V_real(1:window_size)));

array_for_table(1:window_size, 7) = round(cumsum(V_imaginary(1:window_size)));

% Insert magnitude and phase into columns of table (no rounding)

array_for_table(1:window_size, 8) = mag(1:window_size);

array_for_table(1:window_size, 9) = phase_angle_deg(1:window_size);

end

end

for n = 1:length(t_input)-1

delta = phase_angle_deg(n+1) - phase_angle_deg(n);

estimated_freq(n) = (delta / ((2 * 180) / fs));

end

valid = estimated_freq > 40 & estimated_freq < 80;

invalid = ~valid;

if sum(valid) < 2

error('Not enough valid points for interpolation');

end

estimated_freq(invalid) = interp1(find(valid), estimated_freq(valid), find(invalid), 'linear', 'extrap');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot: Input and Phase Angle AFTER FILTER IS APPLIED TO SIGNAL

figure;

subplot(2,1,1);

plot(t_input, estimated_freq, 'b', 'LineWidth', 1);

title('Phasor Magnitude');

xlabel('Time (s)');

xlim([t_input(window_size * 2), t_input(end)]);

if(datapoints)

xlim([t_input(1), t_input(end)]);

xticks(0:T_input:t_input(end))

xtickformat('%.4f'); % shows more precise decimals

end

ylim([min(estimated_freq(window_size * 2:end)) - 2, max(estimated_freq(window_size * 2:end)) + 2])

ylabel('Magnitude');

grid on;
```

