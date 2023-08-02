function [Mean, Sigma, Max]=stdfit4(Tcorr,titulo)

plotit=0;
if nargin==2, plotit=1;end

nbins=200;
%[par,chi2]=hfitg(Tcorr,nbins,0,-20,20);
%Mean=par(1);Sigma=par(2);Max=par(3);

%Limit=Sigma*sqrt(2*log(Max));
%Limit=Sigma*0.7;
%if Sigma<7
%[par,chi2]=hfitg(Tcorr,nbins,0,-Limit+Mean,Limit+Mean);
%[par,chi2]=hfitg(Tcorr,nbins,0,-1.0,1.0);
%Mean=par(1);Sigma=par(2);Max=par(3);
%Limit=Sigma*0.7;
%end

Bin=0.1; % cm
[N,X]=histf(Tcorr,-1.0:Bin:1.0);
[par,chi2]=hfitg(X,N,plotit);
%[par,chi2]=hfitg(Tcorr,nbins,plotit,-Limit+Mean,Limit+Mean);
Mean=par(1);Sigma=par(2);Max=par(3);

if plotit 
	yaxis(1,par(3)*2); logy
	xaxis(-3,3);
	title(titulo);
	drawnow
end

return