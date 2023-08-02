function Y = MERT_Q_for_fit(par, X)

global alphaD alphaQ

A1 = 0; D = 0; F = 0; H = 0; A2 = 0; G = 0;

A  = par(1);
if(length(par)>=2), A1 = par(2); end
if(length(par)>=3), D  = par(3); end
if(length(par)>=4), F  = par(4); end
if(length(par)>=5), H  = par(5); end
if(length(par)>=6), A2 = par(6); end
if(length(par)>=7), G  = par(7); end

Y = 10^20 * Q_MERT(X, alphaD, alphaQ, A, A1, D, F, H, A2, G);  %A^2

end

