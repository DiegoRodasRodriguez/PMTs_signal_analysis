function MinFun = fitter_Exp2Conv(par)

global Ndata Tdata sigmaNdata;

Nmodel = spectrum_Exp2Conv(par, Tdata);

MinFun  = (Nmodel-Ndata)./sigmaNdata;

end