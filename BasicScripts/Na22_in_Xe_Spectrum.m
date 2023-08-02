function Y = Na22_in_Xe_Spectrum(par, X, title)

Amp511    = par(1);
Mean511   = par(2);
Sigma511  = par(3);
AmpBck    = par(4);
MeanBck   = par(5);
SigmaBck  = par(6);
AmpStep   = par(7);
ThStep    = par(8);
SlopeStep = par(9);

Step      = AmpStep*(1-1./(1+exp(-SlopeStep*(X-ThStep))));
G511      = Amp511*exp(-(X-Mean511).^2/(2*Sigma511^2));
GBck      = AmpBck*exp(-(X-Mean511+30).^2/(2*Sigma511^2));

Y = Step + G511 + GBck;
Y = Y';

if nargin>2    
    plot(X, Step,'r'); lastline('Linewidth',2);
    plot(X, G511,'b'); lastline('Linewidth',2); hold on;
    plot(X, GBck, 'g'); lastline('Linewidth',2);
    plot(X, Y); lastline('Linewidth',2); lastline('color', 'k');
end

end

