function [Tuncorr,corr]=uncorr_D(Tm, Ts, order, titulo)
%function [Tuncorr,p]=uncorr(Tmain, Tsecond,plotit)
% igual ao uncorr mas sem limites
plotit=1;
if nargin==3, plotit=0; end

CoefCorr=corrcoef(Ts,Tm); corr=CoefCorr(1,2);
p=polyfit(Ts,Tm,order);
Tuncorr=Tm-polyval(p,Ts);

if plotit
figure;
subplot(1,2,1);
plot(Ts,Tm,'.b');hold on;plot(Ts,polyval(p,Ts),'ro');
yaxis(mean(Tm)-std(Tm),mean(Tm)+std(Tm));
title(titulo);
subplot(1,2,2);
plot(Ts,Tuncorr,'.k');
yaxis(mean(Tuncorr)-std(Tuncorr),mean(Tuncorr)+std(Tuncorr));
drawnow
end

