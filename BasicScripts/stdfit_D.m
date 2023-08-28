function [Mean, Sigma, Max]=stdfit_D(Tcorr,titulo)

TDCbinSize=50;  %Make global

plotit=0;
if nargin==2, plotit=1;end

% 1 ns window
maxtime=1000;   
mintime=-1000;
nbins=floor((maxtime-mintime)/TDCbinSize);

[par,chi2]=hfitg(Tcorr,nbins,0,mintime,maxtime);
Mean=par(1);Sigma=par(2);Max=par(3);

Limit = Sigma*1.5;  %Ajusta en 3 sigmas

if Sigma<600
    
nbins=floor(2*Limit/TDCbinSize);
[par,chi2]=hfitg(Tcorr,nbins,0,-Limit+Mean,Limit+Mean);
Mean=par(1);Sigma=par(2);Max=par(3);
Limit=Sigma*1.5;
end

Bin=TDCbinSize; % ps
[N,X]=histf(Tcorr,-Limit+Mean:Bin:Limit+Mean);
[par,chi2]=hfitg(X,N,plotit);
Mean=par(1);Sigma=par(2);Max=par(3);

if plotit 
	yaxis(1,par(3)*2); logy
	xaxis(mintime,maxtime);
	title(titulo);
	drawnow
end



return