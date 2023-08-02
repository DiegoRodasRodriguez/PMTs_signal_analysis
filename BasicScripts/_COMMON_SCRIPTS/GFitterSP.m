function MinFun = GFitterSP(par)

global Ndata Qdata sigmaNdata;         %#ok<NUSED> %For fit

Nmodel = SP_spectrum(par, Qdata);

MinFun = Nmodel-Ndata; 
%MinFun  = (Nmodel-Ndata)./sigmaNdata;
end