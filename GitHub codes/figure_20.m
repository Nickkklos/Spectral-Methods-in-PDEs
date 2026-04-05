clear; clc; close all;

%% exact solution and RHS
u_exact = @(x,y) cos(pi*x).*(y.^2 - 1).^2;
f_rhs   = @(x,y) cos(pi*x).*((pi^2+1)*(y.^2 - 1).^2 - 12*y.^2 + 4);

%% mesh parameter
N = 24;   % number of subintervals in each direction

%% build uniform mesh on [-1,1]^2
x1d = linspace(-1,1,N+1);
[X,Y] = meshgrid(x1d,x1d);
p = [X(:), Y(:)];          % node coordinates, size = (N+1)^2 x 2

% triangulation: split each square into 2 triangles
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

numNodes = size(p,1);
numElem  = size(t,1);

%% assemble global matrices
S = sparse(numNodes,numNodes);   % stiffness matrix
M = sparse(numNodes,numNodes);   % mass matrix
b = zeros(numNodes,1);           % load vector

% 3-point quadrature on reference triangle
lambda = [1/6 1/6 2/3;
          1/6 2/3 1/6;
          2/3 1/6 1/6];

for K = 1:numElem
    nodes = t(K,:);
    x = p(nodes,1);
    y = p(nodes,2);

    % area of triangle
    area = abs(det([1 x(1) y(1);
                    1 x(2) y(2);
                    1 x(3) y(3)]))/2;

    % local stiffness matrix for P1 element
    beta  = [y(2)-y(3); y(3)-y(1); y(1)-y(2)];
    gamma = [x(3)-x(2); x(1)-x(3); x(2)-x(1)];
    Se = (beta*beta' + gamma*gamma')/(4*area);

    % local mass matrix
    Me = area/12 * [2 1 1;
                    1 2 1;
                    1 1 2];

    % local load vector
    be = zeros(3,1);
    for q = 1:3
        phi = lambda(q,:)';      % basis values at quadrature point
        xq = phi.'*x;
        yq = phi.'*y;
        fq = f_rhs(xq,yq);
        be = be + (area/3) * fq * phi;
    end

    % assemble
    S(nodes,nodes) = S(nodes,nodes) + Se;
    M(nodes,nodes) = M(nodes,nodes) + Me;
    b(nodes) = b(nodes) + be;
end

%% solve FE system
% weak form: \int grad u · grad v + \int u v = \int f v
A = S + M;
uh = A \ b;

%% exact solution at nodes
ue = u_exact(p(:,1), p(:,2));
err = uh - ue;

% nodal max error
err_inf = norm(err,inf);

% approximate L2 error using mass matrix
err_L2 = sqrt(err' * M * err);

fprintf('N = %d\n', N);
fprintf('nodes per direction = %d\n', N+1);
fprintf('total nodes = %d\n', numNodes);
fprintf('number of triangles = %d\n', numElem);
fprintf('nodal max error = %.6e\n', err_inf);
fprintf('approx L2 error = %.6e\n', err_L2);

%% plots
figure;
trisurf(t, p(:,1), p(:,2), uh, 'EdgeColor', [0.6 0.6 0.6]);
xlabel('x'); ylabel('y'); zlabel('u_h');
title(['P1 FE solution, N = ', num2str(N)]);
view(45,30); axis tight; box on;

figure;
trisurf(t, p(:,1), p(:,2), abs(err), 'EdgeColor', [0.6 0.6 0.6]);
xlabel('x'); ylabel('y'); zlabel('|u_h-u|');
title(['Nodal error, N = ', num2str(N), ...
       ', max err = ', num2str(err_inf,'%.2e')]);
view(45,30); axis tight; box on;