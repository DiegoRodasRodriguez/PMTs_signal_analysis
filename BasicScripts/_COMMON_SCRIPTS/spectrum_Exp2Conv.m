%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%      Fit to the sum of two exponentials convoluted with a gaussian      %
%            (suited for 2nd continuum in Xe above 1bar if not            %
%                   accounting for the formation time)                    %
%                       DGD, 05/01/2021                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = spectrum_Exp2Conv(par, X, tit)

A1     = par(1);
tau1   = par(2);
A3     = par(3);
tau3   = par(4);
sigma0 = par(5);
t0     = par(6);

Y3 = A3/(2*tau3)*exp(0.5*(sigma0/tau3)^2)*exp(-X/tau3).*(1-erf(-(X - t0 + (sigma0^2)/tau3)/(sqrt(2)*sigma0)));
Y1 = A1/(2*tau1)*exp(0.5*(sigma0/tau1)^2)*exp(-X/tau1).*(1-erf(-(X - t0 + (sigma0^2)/tau1)/(sqrt(2)*sigma0)));

Y  = Y1 + Y3;

if (nargin==3)
    plot(X, Y1, 'b'); lastline('Linewidth',1);
    plot(X, Y3, 'r'); lastline('Linewidth',1);
    plot(X, Y,  'k'); lastline('Linewidth',1);
    legend('data', 'component1', 'component2', 'all');
    title(tit); box on;
end

end

