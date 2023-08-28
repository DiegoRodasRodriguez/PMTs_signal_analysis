function MinFun = GFitterMERT(par)

global Ydata Xdata sigmaYdata

Ymodel       = MERT_Q_for_fit(par, Xdata);
Ymodel_      = Ymodel(Xdata<2);
Ydata_       = Ydata(Xdata<2);
sigmaYdata_  = sigmaYdata(Xdata<2);
%MinFun  = (Ymodel_ - Ydata_).^2; 
MinFun  = sqrt(((Ymodel_ - Ydata_).^2)./(sigmaYdata_.^2));
end