
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 720; % Sampling frequency (Hz)

fs_input = fs_input * 2;

T_input = 1 / fs_input; % Sampling period (s)

f0_input = 55; % Signal frequency (Hz)

t_input = 0:T_input:(((fs_input/2)/f0_input)/2)/f0_input; % Time vector (0.1 seconds)

Vm_input = 10; % Amplitude

omega_input = 2 * pi * f0_input; % Angular frequency

datapoints = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters FOR FILTER

%#ok<*UNRCH>

fs = fs_input * 2; % Sampling frequency

f0_high_range_value = 75;

f0_low_range_value = 35;

f0 = f0_high_range_value; % Filter Frequency that we are esimating. Set this to anywhere between 40 to 70Hz. 60Hz by default for mid point.

T = 1 / fs; % Sampling period

omega = 2 * pi * f0; % Discrete angular frequency (radians/sample)

f_range = linspace(0, fs, 1000); % Frequency range to plot full range. This is only relevant later on when assigning z = exp(j omega T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input data points

% x = [714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699,714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699, 714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699];

x = Vm_input * sin(omega_input * t_input + pi/18) + Vm_input * sin(2*omega_input * t_input) + Vm_input * sin(3*omega_input * t_input) + Vm_input * sin(4*omega_input * t_input) + Vm_input * sin(5*omega_input * t_input) + Vm_input * sin(6*omega_input * t_input); % Input waveform

datapoints = false;

t_input = 0:T_input:((length(x)-1) * T_input);

estimated_freq = f0 * ones(1, length(t_input));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LES filter generator

filter_choice = 0; % Default initialization

track_active_filter = []; % For error handeling for user selection of filter

% Decide which filter to ACTIVATE. Manually configured filters.

fundamental_filter = true;

second_harmonic_filter = true;

third_harmonic_filter = true;

fourth_harmonic_filter = true;

fifth_harmonic_filter = true;

sixth_harmonic_filter = true;

dc_filter = false;

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

matrix_A = [];

real_values = []; % Real part of filter. H(z) Vp Cos(Theta)

imaginary_values = []; % Imaginary part of filter. H(z) Vp Sin(Theta)

function generate_matrix_A = generate_filter_matrix(k, window, T, f0, fundamental_filter, second_harmonic_filter, third_harmonic_filter, fourth_harmonic_filter, fifth_harmonic_filter, sixth_harmonic_filter, dc_filter)

% Initialize output matrix

num_features = 2 * (fundamental_filter + second_harmonic_filter + third_harmonic_filter + fourth_harmonic_filter + fifth_harmonic_filter + sixth_harmonic_filter) + 2 * dc_filter;

generate_matrix_A = zeros(window, num_features);

% Initial power offset (centered around 0)

z_power = -1 * (window - 1) / 2;

for n = 1:window

array_of_equations = []; % Reset for each row

if fundamental_filter

array_of_equations(end + 1) = sin(2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(2*pi*f0* T * z_power);

end

if second_harmonic_filter

array_of_equations(end + 1) = sin(2*2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(2*2*pi*f0* T * z_power);

end

if third_harmonic_filter

array_of_equations(end + 1) = sin(3*2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(3*2*pi*f0* T * z_power);

end

if fourth_harmonic_filter

array_of_equations(end + 1) = sin(4*2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(4*2*pi*f0* T * z_power);

end

if fifth_harmonic_filter

array_of_equations(end + 1) = sin(5*2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(5*2*pi*f0* T * z_power);

end

if sixth_harmonic_filter

array_of_equations(end + 1) = sin(6*2*pi*f0* T * z_power);

array_of_equations(end + 1) = cos(6*2*pi*f0* T * z_power);

end

if dc_filter

array_of_equations(end + 1) = 1;

array_of_equations(end + 1) = z_power;

end

generate_matrix_A(n, :) = array_of_equations;

z_power = z_power + 1;

end

end

target_array_location = 1;

%-----------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [estimated_freq, phase_angle_deg, mag] = run_freq_estimation(estimated_freq, x, t_input, f0, fs, window, samples, T, fundamental_filter, second_harmonic_filter, third_harmonic_filter, fourth_harmonic_filter, fifth_harmonic_filter, sixth_harmonic_filter, dc_filter, target_array_location)

window_size = samples;

disp(window_size)

phase_angle_deg = zeros(1, length(t_input));

mag = zeros(1, length(t_input));

x_buffer = zeros(1, window_size);

for n = 2:length(t_input)

tolerance = 0.02;

x_buffer = [x_buffer(2:end), x(n)];

new_test_freq = f0 * ones(1, 4);

for k = 1:3

matrix_A = generate_filter_matrix(n, window, T, new_test_freq(k), fundamental_filter, second_harmonic_filter, third_harmonic_filter, fourth_harmonic_filter, fifth_harmonic_filter, sixth_harmonic_filter, dc_filter);

target_pinv_A = pinv(matrix_A);

real_values = target_pinv_A(target_array_location, :);

imaginary_values = target_pinv_A(target_array_location + 1, :);

V_real = real_values .* x_buffer;

V_imaginary = imaginary_values .* x_buffer;

imaginary_part_Vp_cos_theta = sum(V_imaginary);

real_part_Vp_sin_theta = sum(V_real);

phase_angle_deg(n) = atan2(imaginary_part_Vp_cos_theta, real_part_Vp_sin_theta) * 180 / pi;

mag(n) = sqrt(imaginary_part_Vp_cos_theta^2 + real_part_Vp_sin_theta^2);

unwrapped_phase = unwrap(deg2rad(phase_angle_deg));

delta = rad2deg(unwrapped_phase(n) - unwrapped_phase(n-1));

new_estimate = abs(delta) / ((2 * 180) / fs);

if new_estimate > 40 && new_estimate < 70

if abs(new_estimate - f0) <= tolerance

estimated_freq(n) = new_estimate;

break;

end

f0 = new_estimate;

new_test_freq(k + 1) = new_estimate;

else

estimated_freq(n) = f0;

end

end

end

end

tolerance = 0.02;

min_match_ratio = 0.6;

for f0 = f0_high_range_value:-1:f0_low_range_value

% Run tracking for current f0

[temp_estimate, phase_angle_deg, mag] = run_freq_estimation( estimated_freq, x, t_input, f0, fs, window, samples, T, fundamental_filter, second_harmonic_filter, third_harmonic_filter, fourth_harmonic_filter, fifth_harmonic_filter, sixth_harmonic_filter, dc_filter, target_array_location);

% Round to 2 decimals and find dominant frequency

dominant_freq = mode(round(temp_estimate, 2));

% Count how many are within ±tolerance of dominant_freq

match_count = sum(abs(temp_estimate - dominant_freq) <= tolerance);

match_ratio = match_count / length(temp_estimate);

if match_ratio >= min_match_ratio

fprintf("✓ Match found: f0 = %.2f Hz (%.2f%% match to dominant freq %.2f)\n", f0, match_ratio * 100, dominant_freq);

estimated_freq = temp_estimate; % ✅ final estimate set only if condition met

break;

else

fprintf("× f0 = %.2f rejected (%.2f%% match to dominant freq %.2f)\n", f0, match_ratio * 100, dominant_freq);

% Do not assign estimated_freq

end

end

% Optional fallback if no match found

if isempty(estimated_freq)

fprintf("No frequency in range %d–%d Hz met the match condition.\n", f0_low_range_value, f0_high_range_value);

estimated_freq = f0_low_range_value * ones(1, length(t_input)); % fallback

end

% Optional final f0 result

final_estimated_f0 = mode(round(estimated_freq, 2));

disp("Final estimated f0:");

disp(final_estimated_f0);

% After main loop completes

avg_freq = round(mean(estimated_freq), 2); % Rounded to 2 decimals to match logic

replacement_f0 = f0_high_range_value;

for i = 1:length(estimated_freq)

if replacement_f0 < avg_freq

break; % Stop modifying once we reach or drop below avg

end

estimated_freq(i) = replacement_f0;

replacement_f0 = replacement_f0 - 1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot: Input and Phase Angle AFTER FILTER IS APPLIED TO SIGNAL

figure;

subplot(2,1,1);

plot(t_input, estimated_freq, 'b', 'LineWidth', 1);

title('Phasor Magnitude');

xlabel('Time (s)');

xlim([0, t_input(end)]);

ylabel('Magnitude');

grid on;
```