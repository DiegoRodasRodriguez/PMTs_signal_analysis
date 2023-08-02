function MinFun = GFitterSPped(par)

global Ndata Qdata sigmaNdata;         %#ok<NUSED> %For fit

Nmodel = SP_spectrumPed(par, Qdata);

MinFun = (Nmodel-Ndata).^2; 
%MinFun  = ((Nmodel-Ndata).^2)./(sigmaNdata.^2);
end