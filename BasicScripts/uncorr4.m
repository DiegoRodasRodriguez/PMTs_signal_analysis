function [Tuncorr,p]=uncorr4(Tm, A,plotit)
%function [Tuncorr,p]=uncorr(Tmain, Tsecond,plotit)
% igual ao uncorr mas sem limites

if nargin==2, plotit=0; end

CoefCorr=corrcoef(A,Tm); CoefCorr=CoefCorr(1,2)
p=polyfit(A,Tm,2);
%p=polyfit(Ts,Tm,1);
Tuncorr=Tm-polyval(p,A);
CoefCorr=corrcoef(A,Tuncorr);CoefCorr=CoefCorr(1,2)

if plotit
figure;plot(A,Tm,'.b');hold on;plot(A,polyval(p,A),'ro');
plot(A,Tuncorr,'.k');
title(plotit);drawnow
end

