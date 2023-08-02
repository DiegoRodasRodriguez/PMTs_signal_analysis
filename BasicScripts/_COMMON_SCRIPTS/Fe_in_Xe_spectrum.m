function Y = Fe_in_Xe_spectrum(par, X, tit)

global Ndata Qdata sigmaNdata is35peak is47peak; %#ok<NUSED>

Amp   = par(1);
Mean  = par(2);
Sigma = par(3);

Y     = Amp*exp(-(X-Mean).^2/(2*Sigma^2));

if (nargin==3)
    plot(X, Y,'r'); lastline('Linewidth',1); hold on;
    legend('data', 'fit');
    title(tit); box on;
end



