function Y=Am241_in_Xe_Fit(par)

global Q;

Amp26   = par(1);
Amp30   = par(2);
Amp60   = par(3);
Mean30  = par(4);
Mean60  = par(5);
Sigma30 = par(6);
Sigma60 = par(7);

Mean26  = Mean30-4; %FIXME. Use exact value
Sigma26 = Sigma30*sqrt(26/30);

end