
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters

%#ok<*UNRCH>

fs = 720; % Sampling frequency (Hz)

T = 1 / fs; % Sampling period (s)

f0 = 60; % Signal frequency (Hz)

t = 0:T:((fs/f0)/2)/f0; % Time vector (0.1 seconds)

omega = 2 * pi * f0; % Discrete angular frequency (radians/sample)

f_range = linspace(0, fs, 1000); % Frequency range to plot full range. This is only relevant later on when assigning z = exp(j omega T)

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

dc_filter = true;

% Add to tracker only if true

if fundamental_filter, track_active_filter(end + 1) = 1; end

if second_harmonic_filter, track_active_filter(end + 1) = 2; end

if third_harmonic_filter, track_active_filter(end + 1) = 3; end

if fourth_harmonic_filter, track_active_filter(end + 1) = 4; end

if fifth_harmonic_filter, track_active_filter(end + 1) = 5; end

if dc_filter, track_active_filter(end + 1) = 6; end

%-----------------------------------------------------------

% Decide which filter to TARGET/APPLY (Only activate one)

filters = {"fundamental", "2ndHarmonic", "3rdHarmonic", "4thHarmonic", "5thHarmonic", "DC"};

filter_choice = 1; % FUNDAMENAL

% filter_choice = 2; % 2nd Harmonic

% filter_choice = 3; % 3rd Harmonic

% filter_choice = 4; % 4th Harmonic

% filter_choice = 5; % 5th Harmonic

% filter_choice = 6; % DC

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

if dc_filter, samples = samples + 2; end

window = 7; % In case you want custom window size, can change it here manually

matrix_A = [];

z_power = -1 * (window - 1)/2; % Assumes window will always be an odd number so that the center becomes V0, then even amouns before and after.

real_values = []; % Real part of filter. H(z) Vp Cos(Theta)

imaginary_values = []; % Imaginary part of filter. H(z) Vp Sin(Theta)

array_of_equations = []; % Reset at the beginning of each loop

% t is originally a vector

new_f_over_time = zeros(size(t));

new_f_over_time_imag = zeros(size(t));

for idx = 1:length(t)

t_i = t(idx);

z_power = -1 * (window - 1)/2;

matrix_A = [];

for n = 1:window

array_of_equations = [];

if fundamental_filter

array_of_equations(end + 1) = sin(omega * t_i);

array_of_equations(end + 1) = cos(omega * t_i);

array_of_equations(end + 1) = 2*pi*t_i*cos(omega * t_i);

array_of_equations(end + 1) = -2*pi*t_i*sin(omega * t_i);

array_of_equations(end + 1) = -2*((pi*t_i)^2)*sin(omega * t_i);

array_of_equations(end + 1) = -2*((pi*t_i)^2)*cos(omega * t_i);

end

if dc_filter

array_of_equations(end + 1) = 1;

array_of_equations(end + 1) = 1;

end

matrix_A(n, :) = array_of_equations;

z_power = z_power + 1;

end

A_left_pinv = pinv(matrix_A);

% Find index of target filter

target_array_location = 1;

for i = 1:length(filters)

if strcmp(filters{i}, target_filter)

break;

end

if ismember(i, track_active_filter)

target_array_location = target_array_location + 2;

end

end

delta_f_real = A_left_pinv(4,:) ./ A_left_pinv(1,:); % Assuming this is meaningful

new_f_real = delta_f_real + f0;

delta_f_imag = A_left_pinv(3,:) ./ A_left_pinv(2,:); % Assuming this is meaningful

new_f_imag = delta_f_imag + f0;

disp(mean(new_f_real));

% Let's just take the **mean** new_f as the summary value for time idx

new_f_over_time(idx) = mean(new_f_real);

new_f_over_time_imag(idx) = mean(new_f_imag);

end

% Plotting

figure;

plot(t, new_f_over_time);

xlabel('Time (s)');

ylabel('New f (Hz)');

title('Frequency shift over time');

grid on;

% Plotting

figure;

plot(t, new_f_over_time_imag);

xlabel('Time (s)');

ylabel('New f (Hz)');

title('Frequency shift over time');

grid on;
```