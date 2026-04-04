% smooth function vs rough function
% set up grid and differential matrix
N = 24;
h = 2*pi/N;
x = h*(1:N)';   %Periodic grid: excludes 0, includes 2*pi     

% Fourier 
col = zeros(N,1);
k = (1:N-1)';
col(2:N) = 0.5*(-1).^k .* cot(k*h/2); %this formula assumes even N
D = toeplitz(col, [col(1); col(N:-1:2)]);

%rough periodic function:
v = max(0, 1-abs(x-pi)/2);
clf
subplot(2,2,1), plot(x,v,'.-')
axis([0,2*pi,-0.5,1.5]), grid on, title('function')
subplot(2,2,2), plot(x,D*v,'.-')
axis([0,2*pi,-1,1]), grid on, title('spectral derivative')

%smooth periodic function"
v = exp(sin(x));
vprime = cos(x).*v;
subplot(2,2,3), plot(x,v,'.-')
axis([0,2*pi,0,3]), grid on
subplot(2,2,4), plot(x,D*v,'.-')
axis([0,2*pi,-2,2]), grid on


% Discrete infinity-norm error on the grid
err_spec = norm(D*v - vprime, inf);
subplot(2,2,4);
text(1.2, 1.6, ['max error = ' num2str(err_spec), ' and N is ' num2str(N)]);



%Pseudospectral derivative vs centered finite difference
% Use the same smooth periodic test function so that the comparison is fair.
figure
N2 = 24;
h2 = 2*pi/N2;
x2 = h2*(1:N2)';

col = zeros(N2,1);
k = (1:N2-1)';
col(2:N2) = 0.5*(-1).^k .* cot(k*h2/2); %this formula assumes even N
D2 = toeplitz(col, [col(1); col(N2:-1:2)]);

v2 = exp(sin(x2));
dv_exact = cos(x2) .* v2;

% Fourier pseudospectral derivative
dv_spec = D2 * v2;
err_spec2 = norm(dv_spec - dv_exact, inf);

subplot(2,2,1), plot(x,dv_spec,'.-')
axis([0,2*pi,-2,2]), grid on, title('Pseudospectral method derivative')
text(1.2, 1.6, ['max error = ' num2str(err_spec2), ', and N is ' num2str(N2)]);



M = 24;
h = 2*pi/M;
x = h*(1:M)';
v = exp(sin(x));
vprime = cos(x).*v;

% Second-order centered finite difference
vp_fd2 = (circshift(v,-1) - circshift(v,1)) / (2*h);
err2 = norm(vp_fd2 - vprime, inf);


subplot(2,2,2);plot(x,vp_fd2,'.-')
text(1.2, 1.6, ['max error = ' num2str(err2) ', and N is ' num2str(M)]);
axis([0,2*pi,-2,2]), grid on,title('finite difference derivative')
