function Y = DLC_T(par, X, tit)

rho300 = par(1);
K      = par(2);
a      = par(3);

Y = rho300 * exp(K*((1./X).^a - (1/300)^a));

if (nargin==3)
    plot(X, Y, 'r'); lastline('Linewidth',1); hold on;
    title(tit); box on;
end

end

