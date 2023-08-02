function MinFun = DLC_Fitter(par)

global RHOdata Tdata;        %#ok<NUSED> %For fit

RHOmodel = DLC_T(par, Tdata);

MinFun = (RHOmodel-RHOdata).^2;
end