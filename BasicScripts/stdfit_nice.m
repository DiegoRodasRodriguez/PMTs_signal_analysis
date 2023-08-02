function [Mean, Sigma, Max, Tails, Tails300]=stdfit_nice(Tcorr,titulo)

FontSize=8;

plotit=0;
if nargin==2, plotit=1;end


TDCbinSize=50;  %Make global

% 1 ns window
maxtime=1000;   
mintime=-1000;
nbins=floor((maxtime-mintime)/TDCbinSize);
[par,chi2]=hfitg(Tcorr,nbins,0,mintime,maxtime);
Mean=par(1);Sigma=par(2);Max=par(3);

Limit=Sigma*1.5;

if Sigma<600    
nbins=floor(2*Limit/TDCbinSize);
[par,chi2]=hfitg(Tcorr,nbins,0,-Limit+Mean,Limit+Mean);
Mean=par(1);Sigma=par(2);Max=par(3);
Limit=Sigma*1.5;
end

Bin=TDCbinSize; % ps
[N,X]=histf(Tcorr,-Limit+Mean:Bin:Limit+Mean);
[par,chi2]=hfitg(X,N,plotit*0);
Mean=par(1);Sigma=par(2);Max=par(3);

Tcorr=Tcorr(find(abs(Tcorr<1000))); % cut strange events
TailsRight=length(find(Tcorr>3*Sigma))/length(Tcorr);
TailsLeft=length(find(Tcorr<-3*Sigma))/length(Tcorr);
Tails=TailsLeft+TailsRight;

Tails300=length(find(abs(Tcorr)>300))/length(Tcorr);



% nice time distribution plot
if plotit
x=-1000:Bin:1000-Bin;  % nicer?
%x=-1000+4:4:1000-4;		% more accurate
%x=-1000+4:4:2000;		% more accurate
y=[Tcorr; x'];
[N,X]=histf(y,x);
figure;stairs(X,N); hold on;lastline('color','k')
x2=-1.5*Sigma:5:1.5*Sigma;
plot(x2+4/2,g(x2,Mean,Sigma)*Max,'k')
lastline('linewidth',2);logy;
%plot(-3*Sigma:5:3*Sigma,g(-3*Sigma:5:3*Sigma,Mean,Sigma)*Max*4,'--b')
plot(x+4/2,g(x,Mean,Sigma)*Max,'--k','linewidth',1)
yaxis(1,2*max(N));
set(gca,'Xtick',[-1000 -500 -300 0 300 500 1000])
set(gca,'FontSize',FontSize)
title(titulo);
xlabel('Time difference (ps)','FontSize',FontSize+2)
ylabel('Events/50ps','FontSize',FontSize+2)

text(min(xaxis)+diff(xaxis)/20,max(yaxis)/2,['Sigma=' num2str(Sigma)],'FontSize',FontSize);
text(min(xaxis)+diff(xaxis)/20,max(yaxis)/3,['3-Sigma Tails=' num2str(Tails*100) '%'],'FontSize',FontSize);
text(min(xaxis)+diff(xaxis)/20,max(yaxis)/5,['300ps Tails=' num2str(Tails300*100) '%'],'FontSize',FontSize);
text(min(xaxis)+diff(xaxis)/20,max(yaxis)/10,['Events=' num2str(length(Tcorr))],'FontSize',FontSize);

end

if plotit 
	yaxis(1,par(3)*2); logy
	xaxis(mintime,maxtime);
	title(titulo);
	drawnow
end

return