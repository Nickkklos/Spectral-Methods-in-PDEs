%%
%linear BV ODE
N = 16;
[D,x] = cheb(N);

% Solve u'' = exp(4x), u(-1)=u(1)=0  on Chebyshev grid
D2 = D^2;
D2i = D2(2:N,2:N);              % interior rows/cols
f  = exp(4*x(2:N));             % RHS at interior points
ui = D2i \ f;                   % interior solution
u  = [0; ui; 0];                % enforce Dirichlet BC

% Dense grid on [-1,1]
xx = linspace(-1,1,2001)';

% Barycentric weights for Chebyshev–Lobatto points
w = (-1).^(0:N)'; 
w(1)   = w(1)/2; 
w(end) = w(end)/2;

% Barycentric interpolation p(xx) from data (x,u)
pp = zeros(size(xx));
for k = 1:length(xx)
    diff = xx(k) - x;
    j = find(abs(diff) < 1e-14, 1);
    if ~isempty(j)
        pp(k) = u(j);                 % exactly at a node
    else
        tmp = w./diff;
        pp(k) = sum(tmp.*u)/sum(tmp); % barycentric formula
    end
end

% True solution
truesolution = (exp(4*xx) - sinh(4)*xx - cosh(4))/16;

% Plot
clf
plot(x,u,'.','markersize',14); hold on
plot(xx,pp,'r-','linewidth',1);
plot(xx,truesolution,'k--','linewidth',1);
grid on
legend('grid values','barycentric interpolant','true solution','location','best')
title(['Chebyshev Differential Matrix with N = ' num2str(N)])
error=norm(truesolution-pp,inf)
text(-0.8,-2,['max error is ',num2str(error)],'FontSize', 19)
xlabel x, ylabel u