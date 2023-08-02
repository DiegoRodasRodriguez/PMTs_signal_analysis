function MinFun=gChi2(par)

global Xdata Ydata sigmaData;

Ymodel = par(3)*exp(-(Xdata-par(1)).^2/(2*par(2)^2));

MinFun  = ((Ymodel-Ydata).^2)./(sigmaData.^2);

end