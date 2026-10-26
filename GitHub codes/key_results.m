function key_results
%KEY_RESULTS Print the numbers quoted under "Key results" in the README.
%   Run from this folder (it calls cheb.m). It takes a few seconds and
%   makes no plots. The three computations repeat figure_19.m,
%   figure_20.m and figure_21.m with the plotting removed.

N = 24;
[errCheb, nCheb] = neumann_chebyshev_error(N);
[errP1, nP1] = neumann_p1_error(N);

fprintf('\n-Lap(u) + u = f on [-1,1]^2, homogeneous Neumann conditions,\n');
fprintf('exact solution u = cos(pi*x)*(y^2-1)^2 (figures 19 and 20)\n');
fprintf('  Chebyshev collocation, N = %d (%d nodes): max nodal error %.2e\n', ...
        N, nCheb, errCheb);
fprintf('  P1 finite elements, %d x %d squares (%d nodes): max nodal error %.2e\n', ...
        N, N, nP1, errP1);

Nlist = 25:2:39;
errs24 = disk_poisson_errors(24, Nlist);
errs18 = disk_poisson_errors(18, Nlist);

fprintf('\nLap(u) = f on the unit disk, Fourier-Chebyshev (figure 21)\n');
fprintf('  max error against the Chebyshev order N\n');
fprintf('     N      M = 24      M = 18\n');
for k = 1:numel(Nlist)
    fprintf('  %4d    %.2e    %.2e\n', Nlist(k), errs24(k), errs18(k));
end
fprintf('\n');
end


function [err, nNodes] = neumann_chebyshev_error(N)
% Chebyshev collocation for -Lap(u) + u = f with homogeneous Neumann
% conditions, as in figure_19.m.

[D,x] = cheb(N); y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
f = cos(pi.*xx).*((pi^2+1)*(yy.^2-1).^2+4-12*yy.^2);
D2 = D^2;  I = eye(N+1);
L = -kron(I,D2) - kron(D2,I) + kron(I,I);

Dx = kron(D,I);
Dy = kron(I,D);

left   = 1:(N+1);
right  = N*(N+1)+(1:(N+1));
top    = 1:(N+1):(N+1)^2;
bottom = (N+1):(N+1):(N+1)^2;

% u_x = 0 on the two x-boundaries (including the corners),
% u_y = 0 on the two y-boundaries
L(left,:)  = Dx(left,:);
L(right,:) = Dx(right,:);
L(top(2:N),:)    = Dy(top(2:N),:);
L(bottom(2:N),:) = Dy(bottom(2:N),:);

bc = [left, right, top(2:N), bottom(2:N)];
f(bc) = 0;

u = L\f;
u_exact = cos(pi.*xx).*(yy.^2-1).^2;

err = max(abs(u - u_exact));
nNodes = numel(u);
end


function [err, nNodes] = neumann_p1_error(N)
% P1 finite elements on a uniform N-by-N mesh of squares, each split into
% two triangles, for the same problem, as in figure_20.m.

u_exact = @(x,y) cos(pi*x).*(y.^2 - 1).^2;
f_rhs   = @(x,y) cos(pi*x).*((pi^2+1)*(y.^2 - 1).^2 - 12*y.^2 + 4);

x1d = linspace(-1,1,N+1);
[X,Y] = meshgrid(x1d,x1d);
p = [X(:), Y(:)];

t = zeros(2*N*N,3);
elem = 0;
for j = 1:N
    for i = 1:N
        n1 = (j-1)*(N+1) + i;
        n2 = n1 + 1;
        n3 = j*(N+1) + i;
        n4 = n3 + 1;

        elem = elem + 1;
        t(elem,:) = [n1, n2, n4];

        elem = elem + 1;
        t(elem,:) = [n1, n4, n3];
    end
end

nNodes = size(p,1);
S = sparse(nNodes,nNodes);   % stiffness matrix
M = sparse(nNodes,nNodes);   % mass matrix
b = zeros(nNodes,1);         % load vector

% 3-point quadrature on the reference triangle
lambda = [1/6 1/6 2/3;
          1/6 2/3 1/6;
          2/3 1/6 1/6];

for K = 1:size(t,1)
    nodes = t(K,:);
    x = p(nodes,1);
    y = p(nodes,2);

    area = abs(det([1 x(1) y(1);
                    1 x(2) y(2);
                    1 x(3) y(3)]))/2;

    beta  = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    gamma = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    Se = (beta*beta' + gamma*gamma')/(4*area);

    Me = area/12 * [2 1 1;
                    1 2 1;
                    1 1 2];

    be = zeros(3,1);
    for q = 1:3
        phi = lambda(q,:)';
        xq = phi.'*x;
        yq = phi.'*y;
        be = be + (area/3) * f_rhs(xq,yq) * phi;
    end

    S(nodes,nodes) = S(nodes,nodes) + Se;
    M(nodes,nodes) = M(nodes,nodes) + Me;
    b(nodes) = b(nodes) + be;
end

uh = (S + M) \ b;
err = norm(uh - u_exact(p(:,1), p(:,2)), inf);
end


function errs = disk_poisson_errors(M, Nlist)
% Fourier-Chebyshev solver for Lap(u) = f on the unit disk, as in
% figure_21.m. M is the (even) Fourier order; Nlist holds odd Chebyshev
% orders.

dt = 2*pi/M;
t  = dt*(0:M-1)';
Mhalf = M/2;

c = [-pi^2/(3*dt^2) - 1/6, ...
     0.5*(-1).^(2:M) ./ sin(dt*(1:M-1)/2).^2];
Dtt = toeplitz(c);

Z = zeros(Mhalf);
I = eye(Mhalf);
J = [Z I; I Z];

errs = zeros(size(Nlist));

for k = 1:numel(Nlist)
    N = Nlist(k);
    [D, r] = cheb(N);
    Nhalf = (N-1)/2;

    Dsq  = D^2;
    D11  = Dsq(2:Nhalf+1, 2:Nhalf+1);
    D12  = Dsq(2:Nhalf+1, N:-1:Nhalf+2);
    E11  = D(2:Nhalf+1, 2:Nhalf+1);
    E12  = D(2:Nhalf+1, N:-1:Nhalf+2);

    rint = r(2:Nhalf+1);
    R1   = diag(1./rint);
    R2   = diag(1./(rint.^2));

    L = kron(D11 + R1*E11, eye(M)) ...
      + kron(D12 + R1*E12, J) ...
      + kron(R2, Dtt);

    [rr, tt] = meshgrid(rint, t);
    A   = exp(-rr.^5 .* cos(tt));
    B   = exp(-rr.^3 .* cos(tt));
    uex = A - B;

    f = A .* (25*rr.^8 .* cos(tt).^2 ...
            - 24*rr.^3 .* cos(tt) ...
            + rr.^8 .* sin(tt).^2) ...
      + B .* (-9*rr.^4 .* cos(tt).^2 ...
            + 8*rr .* cos(tt) ...
            - rr.^4 .* sin(tt).^2);

    u = L \ f(:);
    errs(k) = max(abs(u - uex(:)));
end
end
