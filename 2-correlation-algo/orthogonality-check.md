
```
% -------- INPUT SECTION --------

% Example with functions:

x = @(t) sin(2*pi*60*t);

y = @(t) 3.5*sin(2*pi*60*t);

domain = [0, 0.1];

% Example with data vectors (comment above if using below)

% x = [1, 2, 3, 4, 5];

% y = [5, 4, 3, -2, -1];

% n = 0:length(x)-1;

% -------- PROCESSING --------

isFunction = isa(x, 'function_handle') && isa(y, 'function_handle');

figure;

if isFunction

% Continuous functions

a = domain(1); b = domain(2);

t = linspace(a, b, 1000);

x_vals = x(t);

y_vals = y(t);

product = x_vals .* y_vals;

inner_product = integral(@(t) x(t).*y(t), a, b);

% Plot each signal

subplot(4,1,1); plot(t, x_vals, 'b', 'LineWidth', 2); grid on;

title('x(t)'); xlabel('t'); ylabel('x');

subplot(4,1,2); plot(t, y_vals, 'r', 'LineWidth', 2); grid on;

title('y(t)'); xlabel('t'); ylabel('y');

subplot(4,1,3);

plot(t, product, 'k--', 'LineWidth', 2); grid on;

hold on;

yline(0, 'm:'); % add a zero reference line

title('x(t)·y(t)'); xlabel('t'); ylabel('x(t)y(t)');

% Continuous functions case visualization - subplot 4

subplot(4,1,4); hold on; grid on;

% Compute product values

product = x_vals .* y_vals;

% Separate positive and negative parts

pos_product = product;

pos_product(pos_product < 0) = 0;

neg_product = product;

neg_product(neg_product > 0) = 0;

% Plot positive area in blue

area(t, pos_product, 'FaceColor', [0 0.5 1], 'EdgeColor', 'none');

% Plot negative area in red

area(t, neg_product, 'FaceColor', [1 0 0], 'EdgeColor', 'none');

% Zero reference line

yline(0, 'k:');

% Title and axis labels

title('x(t)·y(t) (Shaded Area)');

xlabel('t');

ylabel('x(t)·y(t)');

% Calculate areas via integration

positive_area = integral(@(t) max(x(t).*y(t), 0), domain(1), domain(2));

negative_area = integral(@(t) min(x(t).*y(t), 0), domain(1), domain(2));

% Annotate positive and negative areas on the plot

text(domain(1) + 0.05*(domain(2)-domain(1)), max(product)*0.9, ...

sprintf('Positive Area = %.4f\nNegative Area = %.4f', positive_area, negative_area), ...

'FontSize', 11, 'Color', 'black', 'BackgroundColor', 'white', 'EdgeColor', 'none');

else

% Discrete vectors

if length(x) ~= length(y)

error("Vectors must be the same length.");

end

if ~exist('n', 'var')

n = 0:length(x)-1;

end

x_vals = x;

y_vals = y;

product = x .* y;

inner_product = sum(product);

subplot(4,1,1); plot(n, x, 'b-o', 'LineWidth', 2); grid on;

title('x[n]'); xlabel('n'); ylabel('x');

subplot(4,1,2); plot(n, y, 'r-o', 'LineWidth', 2); grid on;

title('y[n]'); xlabel('n'); ylabel('y');

subplot(4,1,3); plot(n, product, 'k--o', 'LineWidth', 2); grid on;

title('x[n]·y[n]'); xlabel('n'); ylabel('x[n]y[n]');

subplot(4,1,4); hold on; grid on;

bar(n, product, 'FaceColor', 'flat'); % Bar chart for area visualization

colormap([1 0 0; 0 0.5 1]); % red for negative, blue for positive

% Color bars manually based on sign

for i = 1:length(product)

if product(i) >= 0

bar(n(i), product(i), 'FaceColor', [0 0.5 1]); % blue

else

bar(n(i), product(i), 'FaceColor', [1 0 0]); % red

end

end

% Zero reference line

yline(0, 'k:');

% Title and labels

title('x[n]·y[n] (Shaded Area)');

xlabel('n');

ylabel('x[n]·y[n]');

% Compute area values

positive_area = sum(product(product > 0));

negative_area = sum(product(product < 0));

% Annotate area totals

text(n(end)+0.5, max(product)*0.9, ...

sprintf('Positive Area = %.2f\nNegative Area = %.2f', ...

positive_area, negative_area), ...

'FontSize', 11, ...

'Color', 'black');

end

% -------- ORTHOGONALITY CHECK --------

tolerance = 1e-10;

if abs(inner_product) < tolerance

status = '✅ x and y are orthogonal.';

else

status = sprintf('❌ x and y are NOT orthogonal.\nInner Product = %.6f', inner_product);

end

% Add unified title and message box at bottom

sgtitle('Orthogonality Check');

% Resize figure to create space for annotation

set(gcf, 'Position', [100, 100, 800, 700]); % [left, bottom, width, height]

% Adjust annotation Y position to be inside view

annotation('textbox', [0.1, 0.01, 0.8, 0.07], ...

'String', status, ...

'HorizontalAlignment', 'center', ...

'VerticalAlignment', 'middle', ...

'FontSize', 12, ...

'FontWeight', 'bold', ...

'EdgeColor', 'none', ...

'Color', 'magenta', ...

'Interpreter', 'none');
```