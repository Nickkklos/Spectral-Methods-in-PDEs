%Differentiation of exp(sin(x))
% Set up periodic grid
N=24;% Must be even number
h=2*pi/N;
x=h*(1:N)';

%smooth periodic function
y=exp(sin(x));yprime=cos(x).*y;

%FFT of y and differentiate in Fourier space
y_hat=fft(y);w_hat=1i*[0:N/2-1,0,-N/2+1:-1]'.*y_hat;
w=real(ifft(w_hat));clf

%rough periodic function
z=max(0,1-abs(x-pi)/2);z_hat=fft(z);

%FFT of z and differentiate in Fourier space
wnew_hat=1i*[0:N/2-1,0,-N/2+1:-1]'.*z_hat;
wnew=real(ifft(wnew_hat));clf

%plot results
subplot(3,2,1), plot(x,z,'.-'),title('function')
axis([0,2*pi,-.5,1.5]), grid on

subplot(3,2,2), plot(x,wnew,'.-'),title('spectral derivative in fft')
axis([0,2*pi,-1,1]), grid on

subplot(3,2,3), plot(x,y,'.-')
axis([0,2*pi,0,3]), grid on

subplot(3,2,4), plot(x,w,'.-')
axis([0,2*pi,-2,2]), grid on

% Error for the smooth example
error_inf = norm(w - yprime, inf);
subplot(3,2,4);
text(1.2, 1.6, ['max error = ' num2str(error_inf), ' and N is ' num2str(N)]);
