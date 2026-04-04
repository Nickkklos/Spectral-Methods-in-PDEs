%% Polynomial interpolation of the Runge-type function
% Polynomial degrees to test
NS=[4,14,24,34];
% Keeping a slightly larger interval than [-1,1] also shows mild extrapolation.
xx= -1.00:0.005:1.00;clf
for j=1:4
    for i =1:2
        N = NS(j);

        % Choose interpolation nodes
        if i==1, s='equispaced points';x=-1+2*(0:N)/N;end
        if i==2, s='Chebyshev points'; x=cos((0:N)*pi/N);end
        subplot(4,2,2*(j-1)+i)

        % The target function is
        u=1./(1+16*x.^2);
        uu=1./(1+16*xx.^2);

        p=polyfit(x,u,N);
        pp=polyval(p,xx);

        subplot(4,2,i+2*(j-1)), plot(xx,pp,'r-'), 
        title(['Polynomial interpoltion using ' s])
        hold on, plot(xx,uu,'b--'), grid on

        legend('Polynomial interpolation', 'True function')
        axis([-1.1,1.1,-1,1.5]);
        max_err=norm(uu-pp,inf);
        text(-0.5,-0.5,['max error = ',num2str(max_err),' N = ',num2str(NS(j))])
    end
end
