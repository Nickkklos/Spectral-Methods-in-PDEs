%2D Elliptic spectral method
function [D,x]=cheb(N)
if N==0,D=0;x=1;return,end
x=cos(pi*(0:N)/N)';
c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=repmat(x,1,N+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D'));
end

N = 24; [D,x] = cheb(N); y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:); yy = yy(:);
f = cos(pi.*xx).*((pi^2+1)*(yy.^2-1).^2+4-12*yy.^2);
D2 = D^2;  I = eye(N+1);
L = -kron(I,D2) - kron(D2,I) + kron(I,I); %construct discretised differential operator

%preparation for Neumann boundary condition
Dx=kron(D,I);
Dy=kron(I,D);

left = 1:(N+1); %left boundary 
right = N*(N+1)+(1:(N+1)); %right boundary
top = 1:(N+1):(N+1)^2; % top boundary
bottom = (N+1):(N+1):(N+1)^2; % bottom boundary

%impose Neumann boundary conditions
L(left,:)  = Dx(left,:);
L(right,:) = Dx(right,:);

%To avoid counting twice we force the four corners to take u_x = 0
L(top(2:N),:)    = Dy(top(2:N),:);
L(bottom(2:N),:) = Dy(bottom(2:N),:);

% homogeneous Neumann: RHS = 0
bc = [left, right, top(2:N), bottom(2:N)];
f(bc) = 0;

%solve the equations
u=L\f; 

% reshape to matrix
[XX,YY] = meshgrid(x,y);
U = reshape(u,N+1,N+1);

% exact solution
u_exact = cos(pi.*xx).*(yy.^2-1).^2;
U_exact = reshape(u_exact,N+1,N+1);

% error
figure
surf(XX,YY,abs(U-U_exact))
shading interp
colorbar
xlabel('x')
ylabel('y')
zlabel('|error|')
title(['Absolute error with N=' ,num2str(N)])