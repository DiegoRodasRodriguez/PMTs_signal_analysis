%Uncorrelation function based on fit to a second order polynomial

function [Tuncorr,p]=uncorr2(Tm, Ts,plotit, verbose)

if nargin==2 
    plotit=0; 
    verbose=0;
end

CoefCorr=corrcoef(Ts,Tm);      
if(verbose==1) 
    CoefCorr=CoefCorr(1,2) 
end;

p=polyfit(Ts,Tm,2);
Tuncorr=Tm-polyval(p,Ts);
CoefCorr=corrcoef(Ts,Tuncorr); 
if(verbose==1) 
    CoefCorr=CoefCorr(1,2) 
end;

if plotit
figure;plot(Ts,Tm,'.b');hold on;plot(Ts,polyval(p,Ts),'ro');
plot(Ts,Tuncorr,'.k');
title(plotit);drawnow
end

