I=find(Tstart<max((max(Tstart)-300),500) & Tstart>30);
Tstart=Tstart(I)*47; % Discard overflow channels
Arpc=Arpc(I); %Appc=Appc(I);
Trpc=Trpc(I); %Tppc=Tppc(I);
TrpcT1=TrpcT1(I); %TppcT1=TppcT1(I);
Tevent=Tevent(I);Ntriggers=Ntriggers(I);
if exist('AS1')
	AS1=AS1(I);
end


meanTstart=mean(Tstart);
Tstart=Tstart-meanTstart;
tt=Tstart(find(abs(Tstart)<20*47));
[N,X]=hist(tt,(-20:20)*47);
%figure,stairs(X,N); logy;
%SigmaStart=std(tt)/sqrt(2)

par=hfitg(tt,30,0);
Mean=par(1);Sigma=par(2);
Limit=Sigma*sqrt(2*log(par(3)));
[par,chi2]=hfitg(tt,30,1,-2*Limit+Mean,2*Limit+Mean); logy
yaxis(1,par(3)*2);

SigmaStart=par(2)/sqrt(2)
meanTstart=par(1);
Tstart=Tstart-meanTstart;

text(min(X),max(N),['sigma/sqrt(2)= ' num2str(SigmaStart) ' ps']);
title('start');
xlabel('ps');ylabel('counts');
xaxis(-1000,1000);
drawnow
