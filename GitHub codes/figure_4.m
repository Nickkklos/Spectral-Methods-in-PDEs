N=128;h=2*pi/N;x=h*(1:N);t=0;dt=h/4;
c=0.2+sin(x-1).^2;
v=exp(-100*(x-1).^2);
x = x(:);          % N×1
c = c(:);          % N×1
v = v(:);          % N×1

col = zeros(N,1);
m = (1:N-1)';                         % offsets
col(2:N) = 0.5*(-1).^m .* cot(m*h/2); % even N formula 
D = toeplitz(col, [col(1); col(N:-1:2)]);



%Output control
tmax = 8; % Target final time
tplot = 0.15; % Time between stored frames

clf; drawnow
plotgap = round(tplot/dt);
dt = tplot/plotgap;
nplots = floor(tmax/tplot);

% --- Backward Euler system matrix A = I + dt*diag(c)*D ---
data = zeros(nplots+1, N);
data(1,:) = v.';     % store initial data
tdata = zeros(nplots+1,1);
tdata(1) = t;

A = eye(N) + dt * (diag(c) * D);
[L,U,P] = lu(A);

%Time stepping
for i =1:nplots
    for n=1:plotgap
        t=t+dt;
        v = U \ (L \ (P*v));
    end

    data(i+1,:) = v.';    
    tdata(i+1) = t;
end
waterfall(x, tdata, data);
axis([0,2*pi,0,tmax,0,3.5])
xlabel('x'); ylabel('t'); zlabel('u');
title('Backward Euler + Fourier differentiation matrix');
view(-70,20); grid off
