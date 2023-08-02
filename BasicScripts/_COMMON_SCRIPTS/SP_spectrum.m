function Y = SP_spectrum(par, X, tit)

nmeanPh  = par(1);
ampQ1    = par(2);
sigmaQ1  = par(3);
meanQ1   = par(4);
sigmaPed = par(5);
meanPed  = par(6);

%assume 4 peaks by default (term exp(-nph) reabsorbed in ampQ1)

G0 = ampQ1/nmeanPh              * exp((-(X - meanPed).^2/(2*sigmaPed^2)));
G1 = ampQ1                      * exp(-(X - (1*meanQ1 + meanPed)).^2/(2*(1*sigmaQ1^2 + sigmaPed^2)));
G2 = ampQ1*1/2* nmeanPh         * exp(-(X - (2*meanQ1 + meanPed)).^2/(2*(2*sigmaQ1^2 + sigmaPed^2)));
G3 = ampQ1*1/6* nmeanPh^2       * exp(-(X - (3*meanQ1 + meanPed)).^2/(2*(3*sigmaQ1^2 + sigmaPed^2)));
G4 = ampQ1*1/24*nmeanPh^3       * exp(-(X - (4*meanQ1 + meanPed)).^2/(2*(4*sigmaQ1^2 + sigmaPed^2)));

G5 = ampQ1*1/120* nmeanPh^4     * exp(-(X - (5*meanQ1 + meanPed)).^2/(2*(5*sigmaQ1^2 + sigmaPed^2)));
G6 = ampQ1*1/720* nmeanPh^5     * exp(-(X - (6*meanQ1 + meanPed)).^2/(2*(6*sigmaQ1^2 + sigmaPed^2)));
G7 = ampQ1*1/5040*nmeanPh^6     * exp(-(X - (7*meanQ1 + meanPed)).^2/(2*(7*sigmaQ1^2 + sigmaPed^2)));

Y = G0 + G1 + G2 + G3 + G4 + G5 + G6 + G7;

if (nargin==3)
    plot(X, G0,'k'); lastline('Linewidth',1); hold on;
    plot(X, G1,'g'); lastline('Linewidth',1);
    plot(X, G2,'r'); lastline('Linewidth',1);
    plot(X, G3,'b'); lastline('Linewidth',1);
    plot(X, G4,'c'); lastline('Linewidth',1);    
    plot(X, G5,'g'); lastline('Linewidth',1);
    plot(X, G6,'r'); lastline('Linewidth',1);
    plot(X, G7,'b'); lastline('Linewidth',1);
    plot(X, Y, 'k'); lastline('Linewidth',2);
    legend('data','pedestal','1ph','2ph','3ph','4ph','5ph','6ph','7ph','all');
    title(tit); box on;
end

end

