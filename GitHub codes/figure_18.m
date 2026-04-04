%% Spectral collocation for a 1D Poisson problem on [-1,1]


clear; clc; close all;

%Dense grid for plotting and error measurement
xx = -1:0.001:1;

% Exact solution satisfying u'' = exp(4x), u(1)=0, u'(-1)=0
uexact = (exp(4*xx) - 4*exp(-4)*(xx - 1) - exp(4)) / 16;

%Solve once for a representative value of N
N = 16;

[D, x] = cheb(N);      % Chebyshev differentiation matrix and nodes
D2 = D^2;              % Second derivative matrix

% Replace the last equation by the Neumann boundary condition at x = -1:
%     u'(-1) = 0
D2(N+1, :) = D(N+1, :);

% Unknowns are u at nodes x(2),...,x(N+1), because u(1)=0 is imposed directly.
A = D2(2:N+1, 2:N+1);

% Right-hand side:
% PDE enforced at x(2),...,x(N)
% Neumann BC enforced in the last row
rhs = [exp(4*x(2:N)); 0];

% Solve for the unknown nodal values and append the Dirichlet value u(1)=0
u = A \ rhs;
u = [0; u];

% Interpolate from Chebyshev nodes to the dense plotting grid
uu = barycheb1(x, u, xx);

%Plot numerical solution and interpolant
clf
subplot(2,1,1)
hold on
plot(x, u, '.', 'MarkerSize', 16)
plot(xx, uu, '-', 'LineWidth', 1.2)
grid on
xlabel('x')
ylabel('u(x)')
title(['Max error = ', num2str(norm(uu - uexact, inf)), ',   N = ', num2str(N)])
legend('Grid values', 'Barycentric interpolant', 'Location', 'best')
xlim([-1 1])

%Convergence study
Ns   = 4:20;                  % Chebyshev polynomial orders
errs = zeros(size(Ns));

for k = 1:length(Ns)
    N = Ns(k);

    [D, x] = cheb(N);
    D2 = D^2;

    % Impose u'(-1)=0 in the last row
    D2(N+1, :) = D(N+1, :);

    % Remove the known Dirichlet node u(1)=0
    A = D2(2:N+1, 2:N+1);
    rhs = [exp(4*x(2:N)); 0];

    u = A \ rhs;
    u = [0; u];

    % Interpolate to the dense grid and compute the max error
    uu = barycheb1(x, u, xx);
    errs(k) = norm(uu - uexact, inf);
end

subplot(2,1,2)
semilogy(Ns, errs, 'o-k', 'LineWidth', 1.5)
grid on
xlabel('N (Chebyshev order)')
ylabel('Max error on dense grid')
title('Spectral convergence for 1D Poisson problem with mixed Dirichlet--Neumann BCs')



