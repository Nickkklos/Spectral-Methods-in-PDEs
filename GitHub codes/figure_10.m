%%
% Accuracy of Chebyshev spectral differentiation for several test functions
Nmax=50; E=zeros(4,Nmax);
for N = 1:Nmax
    [D,x]=cheb(N);

    %first test function
    v=abs(x).^3;vprime=3*x.*abs(x);
    E(1,N)=norm(D*v-vprime,inf);

    %second test function
    v = zeros(size(x));
    vprime = zeros(size(x));
    idx = (x ~= 0);                 % avoid division by zero
    v(idx) = exp(-1./x(idx).^2);
    vprime(idx) = 2*v(idx)./x(idx).^3;
    E(2, N) = norm(D*v - vprime, inf);

    %third test function
    v=1./(1+x.^2);vprime=-2*x.*v.^2;
    E(3,N) = norm(D*v - vprime, inf);

    %fourth test function
    v=x.^10;vprime=10*x.^9;
    E(4,N) = norm(D*v - vprime, inf);
end
titles={'|x|^3','exp(-x^{-2})','1/(1+x^2)','x^{10}'};clf
for iplot=1:4
    subplot(2,2,iplot)
     % Plot error against N on a semi scale
    semilogy(1:Nmax,E(iplot,:),'.-') 
    grid on
    
    axis([1,Nmax,10^-16,10^3]),grid on
    set(gca,'xtick',0:10:Nmax,'ytick',10.^(-15:5:0))
    xlabel N, ylabel error, title(titles(iplot))
end
