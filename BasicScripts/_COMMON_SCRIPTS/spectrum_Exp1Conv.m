%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%           Fit to an exponential convoluted with a Gaussian              %
%                   (suited for 3rd continuum)                            %
%                        DGD, 05/01/2021                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Y = spectrum_Exp1Conv(par, X, tit)

A      = par(1);
tau    = par(2);
sigma0 = par(3);
t0     = par(4);

Y = A/(2*tau)*exp(0.5*(sigma0/tau)^2)*exp(-X/tau).*(1-erf(-(X - t0 + (sigma0^2)/tau)/(sqrt(2)*sigma0)));

if (nargin==3)
    plot(X, Y,'k'); lastline('Linewidth',1);
    legend('data','fit');
    title(tit); box on;
end

end

