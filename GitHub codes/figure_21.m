%% Spectral convergence for Poisson equation on the unit disk

clear; clc; close all;

% Chebyshev orders to test (must be odd)
Nlist = 25:2:39;

% Fourier orders to compare (must be even)
M1 = 24;
M2 = 18;

% Compute errors for each Fourier order
errs1 = compute_disk_poisson_errors(M1, Nlist);
errs2 = compute_disk_poisson_errors(M2, Nlist);

% Plot convergence
figure;

subplot(2,1,1)
semilogy(Nlist, errs1, 'o-', 'LineWidth', 1.6);
grid on
xlabel('Chebyshev order N');
ylabel('Max error');
title(['Spectral convergence for \Delta u = f on the unit disk, M = ' num2str(M1)]);

subplot(2,1,2)
semilogy(Nlist, errs2, 'o-', 'LineWidth', 1.6);
grid on
xlabel('Chebyshev order N');
ylabel('Max error');
title(['Spectral convergence for \Delta u = f on the unit disk, M = ' num2str(M2)]);


% Local function
function errs = compute_disk_poisson_errors(M, Nlist)
% compute_disk_poisson_errors
% Computes the max-norm error for the spectral solution of Delta u = f
% on the unit disk for a fixed Fourier order M and a list of odd
% Chebyshev orders N.

    if mod(M,2) ~= 0
        error('Fourier order M must be even.');
    end

    % Angular grid and Fourier second-derivative matrix
    dt = 2*pi/M;
    t  = dt*(0:M-1)';
    Mhalf = M/2;

    % Fourier spectral second derivative matrix on equispaced grid
    c = [-pi^2/(3*dt^2) - 1/6, ...
         0.5*(-1).^(2:M) ./ sin(dt*(1:M-1)/2).^2];
    Dtt = toeplitz(c);

    % Permutation-type block appearing in the disk method
    Z = zeros(Mhalf);
    I = eye(Mhalf);
    J = [Z I; I Z];

    errs = zeros(size(Nlist));

    for k = 1:numel(Nlist)
        N = Nlist(k);

        if mod(N,2) ~= 1
            error('Chebyshev order N must be odd.');
        end

        % Chebyshev differentiation matrix on [-1,1]
        [D, r] = cheb(N);

        % Number of positive interior radial points retained
        Nhalf = (N-1)/2;

        % First and second radial differentiation blocks
        Dsq  = D^2;
        D11  = Dsq(2:Nhalf+1, 2:Nhalf+1);
        D12  = Dsq(2:Nhalf+1, N:-1:Nhalf+2);

        E11  = D(2:Nhalf+1, 2:Nhalf+1);
        E12  = D(2:Nhalf+1, N:-1:Nhalf+2);

        % Diagonal factors 1/r and 1/r^2
        rint = r(2:Nhalf+1);
        R1   = diag(1./rint);
        R2   = diag(1./(rint.^2));

        % Discrete polar Laplacian on the disk:
        %   u_rr + (1/r)u_r + (1/r^2)u_{theta theta}
        L = kron(D11 + R1*E11, eye(M)) ...
          + kron(D12 + R1*E12, J) ...
          + kron(R2, Dtt);

        % Grid in polar coordinates
        [rr, tt] = meshgrid(rint, t);   % size M x Nhalf

        % Exact solution
        A   = exp(-rr.^5 .* cos(tt));
        B   = exp(-rr.^3 .* cos(tt));
        uex = A - B;

        % RHS f = Delta u corresponding to the exact solution above
        f = A .* (25*rr.^8 .* cos(tt).^2 ...
                - 24*rr.^3 .* cos(tt) ...
                + rr.^8 .* sin(tt).^2) ...
          + B .* (-9*rr.^4 .* cos(tt).^2 ...
                + 8*rr .* cos(tt) ...
                - rr.^4 .* sin(tt).^2);

        % Solve the linear system
        u = L \ f(:);

        % Max pointwise error
        errs(k) = max(abs(u - uex(:)));
    end
end
