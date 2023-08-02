function MinFun = fitter_Exp1Conv(par)

global Ndata Tdata sigmaNdata;

Nmodel = spectrum_Exp1Conv(par, Tdata);

MinFun  = (Nmodel-Ndata)./sigmaNdata;

end