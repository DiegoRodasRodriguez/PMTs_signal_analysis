function Y = SP_spectrumPed(par, X, tit)

meanPed  = par(1);
sigmaPed = par(2);
aPed     = par(3);

Y = aPed * exp((-(X - meanPed).^2/(2*sigmaPed^2)));

if (nargin==3)
    plot(X, Y,'k'); lastline('Linewidth',1); hold on;
    legend('data','fit');
    title(tit); box on;
end

end

