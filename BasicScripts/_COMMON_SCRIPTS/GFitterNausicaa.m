function MinFun = GFitterNausicaa(par)

global Ndata Qdata sigmaNdata is35peak is47peak;         %#ok<NUSED> %For fit

Nmodel = Ba_in_Xe_spectrum(par, Qdata);

%MinFun= (Nmodel-Ndata).^2; %weights to be added
%Adjust fit around the peak
Icut        = (Qdata>27 & Qdata<37)|(Qdata>44 & Qdata<61)|(Qdata>78 & Qdata<90);
Nmodel_     = Nmodel(Icut);
Ndata_      = Ndata(Icut);
sigmaNdata_ = sigmaNdata(Icut);

MinFun = (Nmodel_-Ndata_).^2; 
%MinFun  = ((Nmodel_-Ndata_).^2)./(sigmaNdata_.^2);
end