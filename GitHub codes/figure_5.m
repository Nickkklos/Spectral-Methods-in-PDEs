%% Variable-coefficient transport equation on a periodic grid

%solve u_t+c(x)u_x=0

%Grid and initial data 
N=128;
h=2*pi/N;
x=h*(1:N);
t=0;
dt=h/4;

c=0.2+sin(x-1).^2;
v=exp(-100*(x-1).^2);

%output parameters
tmax=8;
tplot=0.15;
%Adjust dt so that exactly 'plotgap' steps fit into each output interval
plotgap = round(tplot / dt);
dt      = tplot / plotgap;
nplots  = round(tmax / tplot);

%Storage for the solution snapshots
data = zeros(nplots+1, N);
data(1,:) = v.';
tdata = zeros(nplots+1, 1);
tdata(1) = 0;

%ensure column vectors display
x = x(:);
c = c(:);
v = v(:);

%construct finite difference differential matrix
e = ones(N,1);
Dfd = spdiags([-e e], [-1 1], N, N) / (2*h);

%Periodic wrap-around entries
Dfd(1,N) = -1/(2*h);
Dfd(N,1) =  1/(2*h);

%Backward Euler 
A = speye(N) + dt * spdiags(c,0,N,N) * Dfd;


[L,U,P,Q,R] = lu(A);


for i = 1:nplots
    for n = 1:plotgap
        t = t + dt;

        % solve A*v_new = v_old via sparse LU factors
        v = Q * (U \ (L \ (P * (R \ v))));
    end

    data(i+1,:) = v.';     
    tdata(i+1) = t;
end
figure;
waterfall(x, tdata, data);
axis([0,2*pi,0,tmax,0,3.5])
xlabel('x'); ylabel('t'); zlabel('u');
title('FD (central) + Backward Euler');
view(-70,20); grid off
