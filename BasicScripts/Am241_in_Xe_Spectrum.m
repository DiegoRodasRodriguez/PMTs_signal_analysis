function Y = Am241_in_Xe_Spectrum(par, X, title)

Amp26   = par(1);
Amp30   = par(2);
Amp60   = par(3);
Mean30  = par(4);
Mean60  = par(5);
Sigma30 = par(6);
Sigma60 = par(7);
Amp34   = par(8);

if(length(par)==15)
    PolBck     = par(9) + par(10)*X + par(11)*X.^2 + par(12)*X.^3;
    GBck       = par(13)*exp(-(X-par(14)).^2/(2*par(15)^2));
    Bck        = PolBck + GBck;
end
Mean26  = Mean30-4;           
Sigma26 = Sigma30*sqrt(26/30);
Mean34  = Mean30+4;           
Sigma34 = Sigma30*sqrt(34/30);

G30     = Amp30*exp(-(X-Mean30).^2/(2*Sigma30^2));
G26     = Amp26*exp(-(X-Mean26).^2/(2*Sigma26^2));
G34     = Amp34*exp(-(X-Mean34).^2/(2*Sigma34^2));
G60     = Amp60*exp(-(X-Mean60).^2/(2*Sigma60^2));

Y = G26 + G30 + G60 + G34;
if(length(par)==15), Y = Y + Bck; end
Y = Y';
if nargin==3
    plot(X, G26,'r'); lastline('Linewidth',2); hold on;
    plot(X, G34,'g'); lastline('Linewidth',2);
    plot(X, G30); lastline('Linewidth',2);
%    plot(X, G60); lastline('Linewidth',2); lastline('color', 'g');    
    plot(X, Y); lastline('Linewidth',2); lastline('color', 'k');
    if(length(par)==15), plot(X, Bck, 'c'); end
end

end

