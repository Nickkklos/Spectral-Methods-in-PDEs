%%Chebyshev spectral differentiation


% Dense grid for plotting
xx=-1.01:0.001:1;
% Function 
uu=exp(xx).*sin(5*xx);
%exact derivative    
trueprime=exp(xx).*(sin(5*xx)+5*cos(5*xx));
clf

%The number of Chebyshev points
NS = [5,10,15,20];

for i = 1:4
    N=NS(i);
    % Construct Chebyshev differentiation matrix and grid
    [D,x]=cheb(N);

    % Spectral derivative at the Chebyshev nodes
    u=exp(x).*sin(5*x);
    prime=D*u;

    %Maximum nodal error in the spectral derivative
    max_err=norm(D*u-exp(x).*(sin(5*x)+5*cos(5*x)),inf);


    subplot(2,2,i)
    plot(x,prime,'.')
    hold on
    
    %Barycentric interpolation of the nodal derivative values
    w = (-1).^(0:N)'; 
    w(1) = w(1)/2; 
    w(end) = w(end)/2;

    pp = zeros(size(xx));
    for k = 1:length(xx)

        diff = xx(k) - x;
        j = find(abs(diff) < 1e-14, 1);
        if ~isempty(j)

            pp(k) = prime(j);
        else
            tmp = w./diff;
            pp(k) = sum(tmp.*prime)/sum(tmp);
        end
    end

    plot(xx,pp,'r-'), title(['N = ' num2str(N)])
    text(-1.3,-8,['max error is ',num2str(max_err)])
    plot(xx,trueprime,'b-'), grid on
    legend('gird value','Chebyshev interpolant', 'True function')
end
