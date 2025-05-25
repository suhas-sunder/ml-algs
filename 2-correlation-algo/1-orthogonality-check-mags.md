
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

% Plot x(t)

subplot(4,1,1);

plot(t_input, x_vals, 'b', 'LineWidth', 2);

grid on;

title('f(t)');

xlabel('t (s)');

ylabel('x');

% Plot y(t)

subplot(4,1,2);

plot(t_input, y_vals, 'r', 'LineWidth', 2);

grid on;

title('g(t)');

xlabel('t (s)');

ylabel('y');

% Plot product x(t)*y(t)

subplot(4,1,3);

plot(t_input, product, 'k--', 'LineWidth', 2);

hold on;

yline(0, 'm:'); % zero reference line

grid on;

title('f(t)·g(t)');

xlabel('t (s)');

ylabel('f(t)·g(t)');

% Plot positive and negative shaded areas of product

subplot(4,1,4);

hold on; grid on;

pos_product = max(product, 0);

neg_product = min(product, 0);

area(t_input, pos_product, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none');

area(t_input, neg_product, 'FaceColor', [1 0 0], 'EdgeColor', 'none');

yline(0, 'k:');

title('f(t)·g(t)');

xlabel('t (s)');

ylabel('f(t)·g(t)');

% Calculate positive and negative areas via integration

positive_area = integral(@(t) max(x(t).*y(t), 0), a, b);

negative_area = integral(@(t) min(x(t).*y(t), 0), a, b);

% Annotate areas on the plot

text(a + 0.05*(b-a), max(product)*0.9, sprintf('Positive Area = %.4f\nNegative Area = %.4f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'white');

else

% Discrete vectors case

if length(x) ~= length(y)

error("Vectors x and y must have the same length.");

end

x_vals = x;

y_vals = y;

product = x .* y;

inner_product = sum(product);

% Plot f[t]

subplot(4,1,1);

stem(n, x_vals, 'b-o', 'LineWidth', 2);

grid on;

title('f[t]');

xlabel('n');

ylabel('x');

% Plot y[t]

subplot(4,1,2);

stem(n, y_vals, 'r-o', 'LineWidth', 2);

grid on;

title('y[t]');

xlabel('n');

ylabel('y');

% Plot product f[t]*g[t]

subplot(4,1,3);

stem(n, product, 'k--o', 'LineWidth', 2);

grid on;

title('f[t]·g[t]');

xlabel('n');

ylabel('f[t]·g[t]');

% Bar plot with coloring for positive/negative product

subplot(4,1,4);

hold on; grid on;

for i = 1:length(product)

if product(i) >= 0

bar(n(i), product(i), 'FaceColor', [0 0.5 1]); % blue

else

bar(n(i), product(i), 'FaceColor', [1 0 0]); % red

end

end

yline(0, 'k:');

title('f[t]·g[t]');

xlabel('n');

ylabel('f[t]·g[t]');

% Calculate positive and negative sums

positive_area = sum(product(product > 0));

negative_area = sum(product(product < 0));

% Annotate area totals

text(n(end)+0.5, max(product)*0.9, sprintf('Positive Area = %.2f\nNegative Area = %.2f', positive_area, negative_area), 'FontSize', 11, 'Color', 'black');

end

% -------- ORTHOGONALITY CHECK --------

tolerance = 1e-10;

if abs(inner_product) < tolerance

status = '✅ f(t) and g(t) are orthogonal.';

else

status = sprintf('❌ f(t) and g(t) are NOT orthogonal.\nInner Product = %.6f', inner_product);

end

% Add unified title and message box at bottom

sgtitle('Orthogonality Check');

% Resize figure to create space for annotation

set(gcf, 'Position', [100, 100, 800, 700]); % [left, bottom, width, height]

% Display orthogonality status below plots

annotation('textbox', [0.1, 0.01, 0.8, 0.07], 'String', status, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'FontWeight', 'bold', 'EdgeColor', 'none', 'Color', 'magenta');
```