
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

col_names = {'Input','Sine','Input * Sin','Real_Part','Cosine','Input * Cos','Imag_Part','Vp', 'Phase_Angle(deg)'};

array_for_table_columns = length(col_names);

array_for_table_rows = samples;

array_for_table = zeros(array_for_table_rows, array_for_table_columns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Works with both data points and formula

x = [714, 2218, 2314, 1233, -99, -1195, -1699, -1029, 714, 2219, 2314, 1233, -99, -1195, -1699];

% x = Vm_input * sin(omega_input * t_input ); % Input waveform

% datapoints = false;

% This resets the length of t_input if x is not an equation (data points instead)

% Insert x into first column of table

array_for_table(1:samples, 1) = x(1:samples);

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

% Allocate arrays to store angles and magnitude values

phase_angle_deg = zeros(1, length(t_input));

mag = zeros(1, length(t_input));

% Apply the 3-sample phasor magnitude and angle estimator

% This is where we actually take 3 SAMPLES and APPLY THE FILTER

window_size = samples; % Window size is same as sample size

real_values = [0, 0.5, 0.866, 1.0, 0.866, 0.5, 0, -0.5, -0.866, -1.0, -0.866, -0.5];

real_values = [real_values(1:window_size)]; % Scale based on sample size

imaginary_values = [1.0, 0.866, 0.5, 0, -0.5, -0.866, -1.0, -0.866, -0.5, 0, 0.5, 0.866];

imaginary_values = [imaginary_values(1:window_size)]; % Scale based on sample size

V_real = zeros(1, window_size);

V_imaginary = zeros(1, window_size);

x_buffer = zeros(1, window_size); % sliding buffer for x

for n = 1:length(t_input)

% Slide buffer: drop oldest, append new sample

x_buffer = [x_buffer(2:end), x(n)]; % This takes all elements of x_buffer except the first one, destructures the array, then we add the result of x(n) at the end

% Apply weights to current buffer

V_real = real_values .* x_buffer;

V_imaginary = imaginary_values .* x_buffer;

imaginary_part_Vp_cos_theta = sum(V_imaginary)/(0.5 * window_size);

real_part_Vp_sin_theta = sum(V_real)/(0.5 * window_size);

phase_angle_deg(n) = atan2(imaginary_part_Vp_cos_theta, real_part_Vp_sin_theta) * 180 / pi;

mag(n) = sqrt(imaginary_part_Vp_cos_theta^2 + real_part_Vp_sin_theta^2);

if(n == samples)

% Insert real and imaginary filter values into columns of table

array_for_table(1:samples, 2) = real_values(1:samples);

array_for_table(1:samples, 5) = imaginary_values(1:samples);

% Insert calculated values into columns for table after calculation is done.

array_for_table(1:samples, 3) = V_real(1:samples);

array_for_table(1:samples, 6) = V_imaginary(1:samples);

array_for_table(1:samples, 4) = cumsum(V_real(1:samples));

array_for_table(1:samples, 7) = cumsum(V_imaginary(1:samples));

% Insert magnitude and phase into columns of table

array_for_table(1:samples, 8) = mag(1:samples);

array_for_table(1:samples, 9) = phase_angle_deg(1:samples);

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot: Input and Phase Angle AFTER FILTER IS APPLIED TO SIGNAL

figure;

subplot(2,1,1);

plot(t_input, mag, 'b', 'LineWidth', 1);

title('Phasor Magnitude');

xlabel('Time (s)');

if(datapoints)

xlim([t_input(1), t_input(end)]);

xticks(0:T_input:t_input(end))

xtickformat('%.4f'); % shows more precise decimals

end

ylabel('Magnitude');

ylim([0, max(mag) + max(mag) * 0.2]);

grid on;

subplot(2,1,2);

plot(t_input, phase_angle_deg, 'r', 'LineWidth', 1);

title('Phasor Phase Angle');

xlabel('Time (s)');

if(datapoints)

xlim([t_input(1), t_input(end)]);

xticks(0:T_input:t_input(end))

xtickformat('%.4f'); % shows more precise decimals

end

ylabel('Angle (degrees)');

ylim([-180 180]);

yticks(-180:60:180); % Set Y-axis ticks at 50-degree intervals

grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOTS CIRCLE GRAPH WITH REAL AND IMAGINARY VALUES OF FILTERED SIGNAL

% We need to account for every sample that needs to be plotted on real/imaginary plane

% Not only every sample, but every sample for multiple cycles for time t_input

samples_per_cycle = round(fs_input / f0_input);

num_cycles = floor((length(t_input)) / samples_per_cycle);

idx_all = 1:(length(t_input));

% Compute phasor coordinates

phasor_real = mag .* cosd(phase_angle_deg);

phasor_imag = mag .* sind(phase_angle_deg);

% Plot

figure; hold on; axis equal;

% Ideal red dashed circle

theta = linspace(0, 2*pi, 300);

plot(round(max(abs(x))) * cos(theta), round(max(abs(x))) * sin(theta), 'r--', 'LineWidth', 1);

% Origin point

plot(0, 0, 'ko', 'MarkerFaceColor', 'none', 'MarkerSize', 4.5);

% All phasor dots (black hollow circles)

plot(phasor_real(idx_all), phasor_imag(idx_all), 'ko', 'MarkerFaceColor', 'none', 'MarkerSize', 4.5);

grid on;

xlabel('Real Axis');

ylabel('Imaginary Axis');

title('Estimated Phasors on Complex Plane');

hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Frequency Analysis of Estimated Phasor using FFT

phasor_complex = mag .* exp(1j * deg2rad(phase_angle_deg)); % Variable phase

phasor_const_phase = mag; % constant phase

N = length(phasor_complex);

half = N/2 + 1;

Y_const = fft(phasor_const_phase);

Y_const_mag = abs(Y_const);

Y_const_mag = Y_const_mag(1:half);

f = (0:half-1) * fs_input / N;

% Plot

figure;

plot(f, Y_const_mag, 'r-', 'LineWidth', 1);

xlabel('Frequency (Hz)');

xlim([f(1), f(end)]);

ylabel('Magnitude');

title('Estimated Phasor FFT');

legend('Variable Phase','Constant Phase');

grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Customize & print table

my_table = array2table(array_for_table, 'VariableNames', col_names);

% Create figure with a reasonable size

fig = figure('Name', sprintf('Values In Window For First %d Samples', samples), 'NumberTitle', 'off', ...

'Position', [100, 100, 900, 300]); % [left bottom width height]

% Convert numeric data to cell array of strings with full precision

formattedData = compose('%.15g', my_table{:,:}); % 15 digits, avoid sci notation

formattedData = reshape(formattedData, size(my_table{:,:}));

% Get numeric matrix from table so that I can remove decimals from column 8 and 9

numericData = my_table{:,:};

[numRows, numCols] = size(numericData);

% Initialize cell array to hold formatted strings

formattedData = strings(numRows, numCols);

% Format each column

for col = 1:numCols

if col == 8 || col == 9

% No decimals for columns 8 and 9

formattedData(:, col) = compose('%d', round(numericData(:, col)));

else

% Full precision for all other columns

formattedData(:, col) = compose('%.15g', numericData(:, col));

end

end

% Convert to cell array

formattedData = cellstr(formattedData);

formattedData = reshape(formattedData, size(numericData));

% Set column widths

baseWidth = 70;

colWidths = repmat({baseWidth * 1.5}, 1, numCols);

% Display table

t = uitable('Parent', fig, ...

'Data', formattedData, ...

'ColumnName', my_table.Properties.VariableNames, ...

'Units', 'normalized', ...

'Position', [0 0 1 1], ...

'RowName', [], ...

'ColumnWidth', colWidths);
```