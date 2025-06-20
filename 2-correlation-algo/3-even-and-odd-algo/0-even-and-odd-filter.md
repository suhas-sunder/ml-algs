
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters

fs = 1440; % Sampling frequency

T = 1 / fs; % Sampling period

f0 = 120;

samples = fs/f0;

cycles = samples/f0;

f = linspace(0, fs, 1000); % Frequency range for plotting

omega = 2 * pi * f; % Discrete angular frequency (radians/sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = exp(1j * omega * T);

% --- H(z) REAL PART Vp COS(Theta) --- The singal itself can be sine or cos

real_values = sin((2 * pi * f0)*(0:T:cycles));

x_sign = sign(real_values); % -1, 0, or 1

% Even and Odd parts of mapped signal

x_sign_rev = fliplr(x_sign);

x_even = 0.5 * (x_sign + x_sign_rev);

x_odd = 0.5 * (x_sign - x_sign_rev);

real_values = x_even + x_odd;

real_values = real_values(2:samples + 1);

real_values_temp = real_values;

disp("real values");

disp(real_values);

real_values = flip(real_values);

% Length of real and imaginary should be the same, so let's pick real for

% sample length.

N = length(real_values);

H_real = zeros(1,length(z));

% --- H(z) IMAGINARY PART Vp SIN(Theta) --- The singal itself can be sine or cos

imaginary_values = cos((2 * pi * f0)*(0:T:cycles));

x_sign = sign(imaginary_values); % -1, 0, or 1

% Even and Odd parts of mapped signal

x_sign_rev = fliplr(x_sign);

x_even = 0.5 * (x_sign + x_sign_rev);

x_odd = 0.5 * (x_sign - x_sign_rev);

imaginary_values = x_even + x_odd;

imaginary_values = imaginary_values(2:samples + 1);

imaginary_values_temp = imaginary_values;

disp("imag values");

disp(imaginary_values_temp);

imaginary_values = flip(imaginary_values);

H_imaginary = zeros(1,length(z));

% Determine H(z) values for real and imaginary

for n = N:-1:1

H_real = H_real + real_values(n) * z.^(-(n-1));

H_imaginary = H_imaginary + imaginary_values(n) * z.^(-(n-1));

% Printing order of H(z) to confirm consecutive order

fprintf("HRe(z) = " + num2str(real_values(n)) + " * z^" + num2str(-(n-1)));

fprintf(" || ");

fprintf("HIm(z) = " + num2str(imaginary_values(n)) + " * z^" + num2str(-(n-1)));

disp(" ");

end

% Calculate factor A

sine_wave_samples = sin(2*pi*(0:N-1)/N); % [0, 0.5, 0.866, 1.0, 0.866, 0.5, 0, -0.5, -0.866, -1.0, -0.866, -0.5];

disp("Sine Wave Samples");

disp(sine_wave_samples);

V_real = zeros(1, N);

V_imaginary = zeros(1, N);

disp(" ")

disp("factor A for X:")

disp(sine_wave_samples .* real_values_temp)

disp("factor A for Y:")

disp(sine_wave_samples .* imaginary_values_temp)

x_buffer = zeros(1, N); % sliding buffer for x

X_For_Factor_A = sum(sine_wave_samples .* real_values);

X_For_Factor_A_temp = sum(sine_wave_samples .* real_values_temp);

disp("Factor A value for X:")

disp(X_For_Factor_A_temp)

Y_For_Factor_A = sum(sine_wave_samples .* imaginary_values);

Y_For_Factor_A_temp = sum(sine_wave_samples .* imaginary_values);

disp("Factor A value for Y:")

disp(Y_For_Factor_A_temp)

factor_A = sqrt(X_For_Factor_A^2 + Y_For_Factor_A^2);

disp("Factor A:")

disp(factor_A)

% Divide real and imaginary values by 1/2 N because the mangnitude is too high.

H_real = H_real/factor_A;

H_imaginary = H_imaginary/factor_A;

% Magnitude of real

mag_H_real = abs(H_real);

% Phase of real

phase_angle_H_real = atan2(imag(H_real), real(H_real)) * 180 / pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude and phase of H

figure;

subplot(2,1,1);

plot(f, mag_H_real, 'b', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Magnitude');

title('Magnitude Response of H(z) Real');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

subplot(2,1,2);

plot(f, phase_angle_H_real, 'r', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Phase (degrees)');

title('Phase Response of H(z) Real');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

yticks(-180:60:180);

ylim([-180 180]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- H2(z) IMAGINARY PART Vp Sin(Theta) ---

mag_H_imaginary = abs(H_imaginary);

phase_angle_H_imaginary = atan2(imag(H_imaginary), real(H_imaginary)) * 180 / pi; % will be zero everywhere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude and phase of H2

figure;

subplot(2,1,1);

plot(f, mag_H_imaginary, 'b', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Magnitude');

title('Magnitude Response of H(z) Imaginary');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

subplot(2,1,2);

plot(f, phase_angle_H_imaginary, 'r', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Phase (degrees)');

title('Phase Response of H(z) Imaginary');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

yticks(-180:60:180);

ylim([-180 180]);
```

![](../images/20250525164350.png)

![](../images/20250525164410.png)