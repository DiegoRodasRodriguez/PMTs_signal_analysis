function MinFun = GFitterNausicaaFe55(par)

global Ndata Qdata sigmaNdata;

Nmodel = Fe_in_Xe_spectrum(par, Qdata);

%MinFun= (Nmodel-Ndata).^2; %weights to be added
%Adjust fit around the peak
Icut        = (Qdata>4 & Qdata<7);
Nmodel_     = Nmodel(Icut);
Ndata_      = Ndata(Icut);
sigmaNdata_ = sigmaNdata(Icut);

MinFun = (Nmodel_ - Ndata_).^2; 
%MinFun  = ((Nmodel_-Ndata_).^2)./(sigmaNdata_.^2);
end