clear all
close all
clc
% Define the energy efficiency function E(x)
E = @(x) ((250 ./ (x + 5)) + 5*(1 + 0.4 * exp(-(x - 3500).^2 / (2 * 500^2))) .* sin(0.05 * x) - 5*(30 ./ (1 + (x - 2500).^2)) + 40 * cos(0.01 * x));

% Define the interval for the search (1500 to 5000 RPM)
x_min = 3000;
x_max = 4500;

% fminbnd finds the minimum, so we negate E(x) to maximize it
[x_opt, E_opt] = fminbnd(@(x) -E(x), x_min, x_max);
%[x_opt, E_opt] = fminunc(@(x) -E(x), x0);

% Display the optimal value and corresponding energy efficiency
disp(['Optimal rotational speed (x): ', num2str(x_opt), ' RPM']);
disp(['Maximum energy efficiency (E(x)): ', num2str(-E_opt), '%']);

% Plot the energy efficiency function to visualize the result
x_vals = linspace(x_min, x_max, 10000);
E_vals = E(x_vals);

figure;
plot(x_vals, E_vals, 'LineWidth', 2);
xlabel('Rotational Speed (RPM)');
ylabel('Energy Efficiency (%)');
title('Energy Efficiency vs Rotational Speed');
grid on;
hold on;
plot(x_opt, -E_opt, 'ro', 'MarkerFaceColor', 'r'); % Highlight the optimum point
legend('Energy Efficiency', 'Optimal Point');