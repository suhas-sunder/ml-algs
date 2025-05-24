
```
clc;

clear;

close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initializing Variables with Parameters

fs = 720; % Sampling frequency

T = 1 / fs; % Sampling period

f = linspace(0, fs, 1000); % Frequency range for plotting

omega = 2 * pi * f; % Discrete angular frequency (radians/sample)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- H(z) REAL PART Vp COS(Theta) ---

z = exp(1j * omega * T);

% Real impulse response values (e.g., 12-point cosine shape)

real_values = [0, 0.5, 0.866, 1.0, 0.866, 0.5, 0, -0.5, -0.866, -1.0, -0.866, -0.5];

imaginary_values = [1.0, 0.866, 0.5, 0, -0.5, -0.866, -1.0, -0.866, -0.5, 0, 0.5, 0.866];

N = length(real_values);

% Calculate H(z) = sum(h[n] * z^(-n))

H_real = zeros(1,length(N));

H_imaginary = zeros(1,length(N));

for n = 1:N

H_real = H_real + real_values(n) * z.^(-(n-1));

H_imaginary = H_imaginary + imaginary_values(n) * z.^(-(n-1));

end

H_real = H_real/(0.5 * N);

H_imaginary = H_imaginary/(0.5 * N);

mag_H1 = abs(H_real);

phase_angle_H1 = atan2(imag(H_real), real(H_real)) * 180 / pi;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude and phase of H

figure;

subplot(2,1,1);

plot(f, mag_H1, 'b', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Magnitude');

title('Magnitude Response of H1(z) Real');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

subplot(2,1,2);

plot(f, phase_angle_H1, 'r', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Phase (degrees)');

title('Phase Response of H1(z) Real');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

yticks(-180:60:180);

ylim([-180 180]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- H2(z) IMAGINARY PART Vp Sin(Theta) ---

mag_H2 = abs(H_imaginary);

phase_angle_H2 = atan2(imag(H_imaginary), real(H_imaginary)) * 180 / pi; % will be zero everywhere

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude and phase of H2

figure;

subplot(2,1,1);

plot(f, mag_H2, 'b', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Magnitude');

title('Magnitude Response of H2(z) Imaginary');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

subplot(2,1,2);

plot(f, phase_angle_H2, 'r', 'LineWidth', 1);

xlabel('Frequency (Hz)');

ylabel('Phase (degrees)');

title('Phase Response of H2(z) Imaginary');

grid on;

xticks(0:60:fs);

xlim([0 fs]);

yticks(-180:60:180);

ylim([-180 180]);
```


![[Pasted image 20250524004117.png]]

![[Pasted image 20250524004137.png]]

