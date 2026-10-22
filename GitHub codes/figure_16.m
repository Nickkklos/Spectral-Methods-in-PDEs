%% 2D Poisson on [-1,1]^2 (Chebyshev collocation) + barycentric plotting
clear; close all; clc
% Uses cheb.m and barycheb1.m from this folder.

N = 14;
[D,x] = cheb(N);  y = x;

% interior grid points
[Xi,Yi] = meshgrid(x(2:N), y(2:N));
xi = Xi(:);  yi = Yi(:);

% Solve Δu = f, with exact u = sin(pi x) sin(pi y)
f = -2*pi^2 * sin(pi*xi) .* sin(pi*yi);

D2 = D^2;
D2 = D2(2:N, 2:N);
I  = eye(N-1);
L  = kron(I, D2) + kron(D2, I);

u = L \ f;

% put into full grid with boundary zeros
uu = zeros(N+1, N+1);
uu(2:N, 2:N) = reshape(u, N-1, N-1);

% error on Chebyshev grid
[XX,YY] = meshgrid(x,y);
u_exact = sin(pi*XX).*sin(pi*YY);
err = norm(uu(:) - u_exact(:), inf);
fprintf('N=%d, max error on Chebyshev grid = %.3e\n', N, err);

%barycentric interpolation to uniform grid for plotting 
xq = (-1:0.04:1)';     % uniform query points
yq = xq;

%interpolate along x for each fixed y-node
Utmp = zeros(N+1, numel(xq));
for j = 1:(N+1)
    Utmp(j,:) = barycheb1(x, uu(j,:).', xq).';
end

% interpolate along y for each fixed xq
Ufine = zeros(numel(yq), numel(xq));
for i = 1:numel(xq)
    Ufine(:,i) = barycheb1(x, Utmp(:,i), yq);
end

% plot
[xxx,yyy] = meshgrid(xq, yq);
figure(1); clf
h = mesh(xxx, yyy, Ufine);

xlabel('x'); ylabel('y'); zlabel('u');
title(sprintf('Chebyshev spectral, N=%d,  max error = %.3e', N, err));
grid on
