
```
% -------- INPUT SECTION --------

clc;

clear;

close all;

%% Initializing Variables with Parameters For INPUT SIGNAL

fs_input = 720; % Sampling frequency (Hz)

T_input = 1 / fs_input; % Sampling period (s)

% Uncomment for continuous-time signals (function handles)

x = @(t) cos(2*pi*60*t);

y = @(t) 3.5*sin(2*pi*60*t);

% Uncomment for discrete signals (data vectors)

% x = [1, 2, 3, 4, 5];

% y = [5, 4, 3, -2, -1];

% Determine signal type

isFunction = isa(x, 'function_handle') && isa(y, 'function_handle');

if isFunction

% Continuous-time signals

% Define a fine time vector over one sample period (adjust as needed)

t_input = 0:T_input:0.1;

else

% Discrete signals

% Define discrete sample indices and corresponding time vector

n = 0:length(x)-1;

t_input = n * T_input;

end

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

text(a * .8, -max(product)*0.9, sprintf('Positive Area = %.4f\nNegative Area = %.4f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'yellow');

else

% Compute values

x_vals = x;

y_vals = y;

product = x .* y;

inner_product = sum(product);

% Combined plot: x[n] and y[n]

subplot(2,1,1);

plot(n, x_vals, 'k-o', 'LineWidth', 2); hold on;

plot(n, y_vals, 'm-o', 'LineWidth', 2);

grid on;

title('Plotting data points');

xlabel('n');

ylabel('Amplitude');

legend('function 1', 'function 2');

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

text(n(end), -max(product)*0.9, sprintf('Positive Area = %.2f\nNegative Area = %.2f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'yellow');

end

% -------- ORTHOGONALITY CHECK --------

tolerance = 1e-10;

if abs(inner_product) < tolerance

status = 'Both functions are orthogonal.';

else

status = sprintf('Both functions are not orthogonal.\nInner Product = %.6f', inner_product);

end

% Add unified title and message box at bottom

sgtitle('');

% Resize figure to create space for annotation

set(gcf, 'Position', [100, 100, 800, 700]); % [left, bottom, width, height]

% Display orthogonality status below plots

annotation('textbox', [0, 0.01, 0.4, 0.95], 'String', status, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', 'magenta');
```


![](../images/20250524224917.png)