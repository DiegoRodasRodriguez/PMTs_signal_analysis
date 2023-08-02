function Y = Ba_in_Xe_spectrum(par, X, tit)

global Ndata Qdata sigmaNdata is35peak is47peak; %#ok<NUSED>

Amp31   = par(1);
Mean31  = par(2);
Sigma31 = par(3);
Amp51   = par(4);   %may include as well Xe K-alpha at 29.8keV (small contribution)
Mean51  = par(5);
Sigma51 = par(6);
Amp81   = par(7);
Mean81  = par(8);
Sigma81 = par(9);

Amp35   = par(10);
Amp47   = par(11);

if(is35peak)
    G35 = Amp35*exp(-(X-Mean31*36/31).^2/(2*(Sigma31^2)*36/31)); %FIXME (36 works better than 35)
else
    G35 = 0;
end

Sigma51 = Sigma31*sqrt(51/31); %FIXME. Fixme Mean 51?

if(is47peak)
    G47 = Amp47*exp(-(X-Mean51*47/51).^2/(2*(Sigma51^2)*47/51));
else
    G47 = 0;
end

G31     = Amp31*exp(-(X-Mean31).^2/(2*Sigma31^2));
G51     = Amp51*exp(-(X-Mean51).^2/(2*Sigma51^2));
G81     = Amp81*exp(-(X-Mean81).^2/(2*Sigma81^2));

Y = G31 + G35 + G47 + G51 + G81;

if (nargin==3 && ~is35peak && ~is47peak)
    plot(X, G31,'r'); lastline('Linewidth',1); hold on;
    plot(X, G51,'g'); lastline('Linewidth',1);
    plot(X, G81,'c'); lastline('Linewidth',1);
    plot(X, Y,  'b'); lastline('Linewidth',1);
    legend('data','31keV','51keV','81keV','sum');
    title(tit); box on;
end
if (nargin==3 && is35peak && ~is47peak)
    plot(X, G31,'r'); lastline('Linewidth',1); hold on;    
    plot(X, G35,'k'); lastline('Linewidth',1);
    plot(X, G51,'g'); lastline('Linewidth',1);
    plot(X, G81,'c'); lastline('Linewidth',1);
    plot(X, Y,  'b'); lastline('Linewidth',1);
    legend('data', '31keV', '35keV', '51keV', '81keV', 'sum');
    title(tit); box on;
end
if (nargin==3 && is35peak && is47peak)
    plot(X, G31,'r'); lastline('Linewidth',1); hold on;    
    plot(X, G35,'k'); lastline('Linewidth',1);
    plot(X, G47,'m'); lastline('Linewidth',1);
    plot(X, G51,'g'); lastline('Linewidth',1);
    plot(X, G81,'c'); lastline('Linewidth',1);
    plot(X, Y,  'b'); lastline('Linewidth',1);
    legend('data', '31keV', '35keV', '47keV', '51keV', '81keV', 'sum');
    title(tit); box on;
end

end

