function [Mean, Sigma, Max]=stdfit(Tcorr,titulo)

plotit=0;
if nargin==2, plotit=1;end

nbins=100;
[par,chi2]=hfitg(Tcorr,nbins,0,-1000,3000);
Mean=par(1);Sigma=par(2);Max=par(3);

%Limit=Sigma*sqrt(2*log(Max));
Limit=Sigma*1.5;
if Sigma<600
[par,chi2]=hfitg(Tcorr,nbins,0,-Limit+Mean,Limit+Mean);
Mean=par(1);Sigma=par(2);Max=par(3);
Limit=Sigma*1.5;
end

Bin=10; % ps
[N,X]=histf(Tcorr,-Limit+Mean:Bin:Limit+Mean);
[par,chi2]=hfitg(X,N,plotit);
%[par,chi2]=hfitg(Tcorr,nbins,plotit,-Limit+Mean,Limit+Mean);
Mean=par(1);Sigma=par(2);Max=par(3);

if plotit 
	yaxis(1,par(3)*2); logy
	xaxis(-1000,3000);
	title(titulo);
	drawnow
end

return