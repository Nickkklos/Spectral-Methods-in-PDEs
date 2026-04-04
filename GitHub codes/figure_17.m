%% 2D Poisson on [-1,1]^2 with homogeneous Dirichlet BC
% Finite Difference (2nd order) version of your spectral code

N = 128;                 % number of subintervals in each direction
h = 2/N;
x = (-1:h:1)';           % (N+1) points
y = x;

% interior grid (exclude boundary)
xi = x(2:N);
yi = y(2:N);
[XXi,YYi] = meshgrid(xi, yi);

% RHS for Δu = f
f = -2*pi^2 * sin(pi*XXi) .* sin(pi*YYi);
f = f(:);                % vectorise, length (N-1)^2

% 1D second-derivative FD matrix on interior points: size (N-1)x(N-1)
e  = ones(N-1,1);
T  = spdiags([e -2*e e], [-1 0 1], N-1, N-1) / h^2;   % approximates d^2/dx^2

I  = speye(N-1);
A  = kron(I, T) + kron(T, I);     % approximates Δ = dxx + dyy

% solve
u = A \ f;

% reshape back to grid (including boundaries)
uu = zeros(N+1, N+1);
uu(2:N, 2:N) = reshape(u, N-1, N-1);

% exact solution on grid
[XX,YY] = meshgrid(x,y);
u_exact = sin(pi*XX) .* sin(pi*YY);

% max error on grid points
err = max(abs(uu(:) - u_exact(:)));

% plot (uniform grid -> interp2 cubic is fine)
[xxx,yyy] = meshgrid(-1:0.04:1, -1:0.04:1);
uuu = interp2(XX,YY,uu,xxx,yyy,"cubic");

mesh(xxx,yyy,uuu),
xlabel x, ylabel y, zlabel u
title(sprintf('FD (2nd order), N = %d, max error = %.3e', N, err))