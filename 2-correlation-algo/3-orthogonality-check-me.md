
```
% -------- INPUT SECTION --------

clc;

clear;

close all;

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 720; % Sampling frequency (Hz)

f0_input = 60; % Signal frequency (Hz)

samples_per_cycle = fs_input/f0_input;

one_cycle = 1/f0_input; % We only want to plot one cycle to check for orthogonality. From 0 to 2Pi/omega = 1/f0

T_input = 1 / fs_input; % Sampling period (s)

% Uncomment for continuous-time signals (function handles)

% x = @(t) cos(2*pi*60*t);

% y = @(t) 3.5*sin(2*pi*60*t);

% Uncomment for discrete signals (data vectors)

x = [1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1 ];

y = [1, 1, 1, -1, -1, -1, -1, -1, -1, 1, 1, 1];

x = x(1:samples_per_cycle);

y = y(1:samples_per_cycle);

t_input = 0:T_input:one_cycle;

n = 0:length(x)-1;

isFunction = isa(x, 'function_handle') && isa(y, 'function_handle');

% -------- PROCESSING --------

figure;

if isFunction

% Continuous functions case

% For inner product integral, integrate over [0, T_input]

a = t_input(1);

b = t_input(end);

% Evaluate signals at t_input points

x_vals = x(t_input);

y_vals = y(t_input);

product = x_vals .* y_vals;

% Compute inner product integral over the interval

inner_product = integral(@(t) x(t).*y(t), a, b);

% Combined plot of x(t) and y(t)

subplot(2,1,1);

plot(t_input, x_vals, 'k', 'LineWidth', 2);

hold on;

plot(t_input, y_vals, 'm', 'LineWidth', 2);

grid on;

title('Plotting both functions');

legend('function 1', 'function 2');

xlabel('t (s)');

ylabel('Amplitude');

% Prepare shaded product components

pos_product = max(product, 0);

neg_product = min(product, 0);

% Combined plot of x(t)*y(t) with shaded areas

subplot(2,1,2);

plot(t_input, product, 'y', 'LineWidth', 2);

hold on;

area(t_input, pos_product, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none');

area(t_input, neg_product, 'FaceColor', [1 0 0], 'EdgeColor', 'none');

yline(0, 'k:');

grid on;

title('Overlapping Area');

xlabel('t (s)');

ylabel('x(t)·y(t)');

% Calculate positive and negative areas via integration

positive_area = integral(@(t) max(x(t).*y(t), 0), a, b);

negative_area = integral(@(t) min(x(t).*y(t), 0), a, b);

% Annotate areas on the plot

text(0, 0, sprintf('Positive Area = %.4f\nNegative Area = %.4f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'yellow');

else

% Compute values

x_vals = x;

y_vals = y;

product = x .* y;

inner_product = sum(product);

% Combined plot: x[n] and y[n]

subplot(4,1,1);

stem(n, x_vals, 'b-o', 'LineWidth', 2);

grid on;

title('Function 1 Data Points');

xlabel('n');

ylabel('x');

% Plot y[t]

subplot(4,1,2);

stem(n, y_vals, 'r-o', 'LineWidth', 2);

grid on;

title('Function 2 Data Points');

xlabel('n');

ylabel('Amplitude');

% Combined plot: x[n]·y[n] with shaded bars

subplot(2,1,2);

hold on; grid on;

for i = 1:length(product)

if product(i) >= 0

bar(n(i), product(i), 'FaceColor', [0 0 0]); % black

else

bar(n(i), product(i), 'FaceColor', [1 0 1]); % magenta

end

end

yline(0, 'k:');

title('Overlapping Area');

xlabel('n');

ylabel('x[n]·y[n]');

% Calculate and annotate area totals

positive_area = sum(product(product > 0));

negative_area = sum(product(product < 0));

text(n(end)-1.5, -max(product)*1.3, sprintf('Positive Area = %.2f\nNegative Area = %.2f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'yellow');

end

% -------- ORTHOGONALITY CHECK --------

tolerance = 1e-10;

if abs(inner_product) < tolerance

status = 'Both functions ARE orthogonal.';

else

status = sprintf('Both functions ARE NOT orthogonal.\nResult = %.6f, which does not = 0!', inner_product);

end

% Add unified title and message box at bottom

sgtitle('Orthogonality Test for 2 Functions (One Cycle)');

% Resize figure to create space for annotation

set(gcf, 'Position', [100, 100, 800, 700]); % [left, bottom, width, height]

% Display orthogonality status below plots

annotation('textbox', [0, 0.01, 0.4, 0.95], 'String', status, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', 'magenta');

fprintf(status);

fprintf(sprintf('\nPositive Area = %.2f\nNegative Area = %.2f\n', positive_area, negative_area));
```


![](../images/20250524224917.png)