% for the June runs

Tcorr=TCorr1;
Arpc=ARpc1;

disp(' ------- apply cuts ---------------------');

%disp('		Tstart cuts')
%MTs=mean(Tstart);[Pars]=hfitg(Tstart,50,0,MTs-250,MTs+250); MTs=Pars(1);SigmaTs=Pars(2);
%[Pars]=hfitg(Tstart,50,1,MTs-0,MTs+50); MTs=Pars(1);SigmaTs=Pars(2);
%	abs(Tstart-MTs)<SigmaTs&... 

disp('		Tracking position cuts')
MPosX1=mean(PosX1);[Pars]=hfitg(PosX1,50,0,MPosX1-10,MPosX1+10); MPosX1=Pars(1);SigmaPosX1=Pars(2);
%abs(PosX1-MPosX1)<SigmaPosX1*1000&... 
MPosY1=mean(PosY1);[Pars]=hfitg(PosY1,50,0,MPosY1-10,MPosY1+10); MPosY1=Pars(1);SigmaPosY1=Pars(2);
%abs(PosY1-MPosY1)<SigmaPosY1*1000&... 

%disp('		Scintilator cuts')
%Aap=(Ch13+Ch14)/2; % comment to run after CAT
MAap=mean(Aap);[Pars]=hfitg(Aap,50,0,MAap-100,MAap); MAap=Pars(1);SigmaAap=Pars(2);

%disp('		Tracking hit strips cut')
%	abs(HitsX1(:,1)-3)<3&...
%	abs(HitsY1(:,1)-3)<3&...
%	abs(HitsX2(:,1)-3)<3&...
%	abs(HitsY2(:,1)-3)<3&...
%	abs(HitsX3(:,1)-3)<3&...
%	abs(HitsY3(:,1)-3)<3&...

[Tc,A]=cut(find(...
	Arpc>100 & Arpc<1000 &...
	abs(Aap-MAap)<SigmaAap*2&...
	Tevent<100 &... 
	abs((HitsX1)-3)<3&...
	abs((HitsY1)-3)<3&...
	abs(PosX1-MPosX1)<SigmaPosX1*2&... 
	abs(PosY1-MPosY1)<SigmaPosY1*2&... 
	TCorr2<5&...
	TCorr3<5&...
	1),Tcorr,Arpc);

%	abs(HitsX2-3)<3&...
%	abs(HitsY2-3)<3&...

%	Tevent>100 & Tevent<200 &...
%	abs(PosX1-MPosX1)<SigmaPosX1*1&... 
%	abs(PosY1-MPosY1)<SigmaPosY1*1&... 

Events=length(Tc)

Tc=Tc+randn(Events,1)*00;

%hfitg(Tcorr,50,0,-500,1000);%logy;yaxis(1,10000);
par=hfitg(Tc,50,0,-500,500);logy;yaxis(1,10000);
Mean=par(1);Sigma=par(2),Max=par(3);
Tc=Tc-Mean;
par=hfitg(Tc,50,1,-3*Sigma,3*Sigma);logy;yaxis(1,10000);
Mean=par(1);Sigma=par(2),Max=par(3);

[N,X]=histf(Tc,-500:20:500);
figure;stairs(X,N/max(N)); hold on
plot(X+10,g(X,Mean,Sigma),'-b');
yaxis(0.011,1);logy;grid on; hold off

[N,X]=histf(Tc,-500:500);
Ns=cumsum(N);
figure; logy; hold on
Gs=Max*cumsum(g(X,Mean,Sigma));Gs=Gs/max(Gs)*max(Ns);
semilogy(X+10,1-Gs/max(Gs),'b');hold on
Gs=Max*cumsum(g(X,Mean,2*Sigma));Gs=Gs/max(Gs)*max(Ns);
semilogy(X+10,1-Gs/max(Gs),'b');hold on
semilogy(X,1-Ns/max(Ns)+1e-7);hold on;
yaxis(0.0027/2,1);grid on; hold off

return

