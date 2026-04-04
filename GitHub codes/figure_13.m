%% Linear two-point BVP solved by finite differences
clear; clc;

% Set up the grid
N = 128;                        % number of subintervals
h = 2 / N;                      % mesh size
x = linspace(-1, 1, N+1)';      % grid points (column vector)

xi = x(2:N);

%Define the right-hand side 
f = exp(4 * xi);

%Build the second-derivative finite difference matrix
e = ones(N-1, 1);
D2fd = spdiags([e, -2*e, e], [-1, 0, 1], N-1, N-1) / h^2;

%Solve the linear system
uiFD = D2fd \ f;

% Append the boundary values to recover the full solution vector
uFD = [0; uiFD; 0];

%  Exact solution

xx = linspace(-1, 1, 2001)';
utrue_xx = (exp(4*xx) - sinh(4)*xx - cosh(4)) / 16;
utrue_x  = (exp(4*x)  - sinh(4)*x  - cosh(4)) / 16;

%Compute the max-norm error on the grid
err = norm(uFD - utrue_x, inf);

%Plot numerical and exact solutions
clf
plot(x, uFD, '.-', 'LineWidth', 1.0, 'MarkerSize', 12); hold on
plot(xx, utrue_xx, 'LineWidth', 1.5)

grid on
xlabel('x')
ylabel('u(x)')
title(['Finite Difference Solution for u'''' = e^{4x},  N = ', num2str(N)])
legend('FD solution on grid', 'Exact solution', 'Location', 'best')

% Put the error text near the lower-left corner of the plot
yl = ylim;
xl = xlim;
text(xl(1) + 0.05*(xl(2)-xl(1)), yl(1) + 0.08*(yl(2)-yl(1)), ...
    ['max error = ', num2str(err, '%.3e')], 'FontSize', 12);