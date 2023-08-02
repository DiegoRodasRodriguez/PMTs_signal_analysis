% for the July runs 2108 to 2122

if 1

clear all
close all

load data\corrall2.mat

TsM=TsM*50;

disp(' ------- tracking ---------------------');


z=200; % position of my chamber
X1Off=16.3; X2Off=14.3; X3Off=13.9;
Y1Off=12.8; Y2Off=12.3; Y3Off=12.4;
Z1=37.5; Z2=314; Z3=426;

PX1_=PX1-X1Off;PX2_=PX2-X2Off;%PX3=PX3-X3Off;
PY1_=PY1-Y1Off;PY2_=PY2-Y2Off;%PY3=PY3-Y3Off;

%figure;subplot(2,3,1);hfitg(PX1,100,2),subplot(2,3,2);hfitg(PX2,100,2),subplot(2,3,3);hfitg(PX3,100,2),
%subplot(2,3,4);hfitg(PY1,100,2),subplot(2,3,5);hfitg(PY2,100,2),subplot(2,3,6);hfitg(PY3,100,2),drawnow

k=(z-Z1)/(Z2-Z1);
x=PX1_+(PX2_-PX1_)*k;
y=PY1_+(PY2_-PY1_)*k;

x=x*4;y=y*4; % convert to mm
%figure;plot(x,y,'.r');axis([-50 50 -50 50]); 

disp(' ------- apply cuts ---------------------');

%disp('		Tstart cuts')
%MTs=mean(Tstart);[Pars]=hfitg(Tstart,50,0,MTs-250,MTs+250); MTs=Pars(1);SigmaTs=Pars(2);
%[Pars]=hfitg(Tstart,50,1,MTs-0,MTs+50); MTs=Pars(1);SigmaTs=Pars(2);
%	abs(Tstart-MTs)<SigmaTs&... 

disp('		Tracking position cuts')
%MPosX1=mean(PX1);[Pars]=hfitg(PX1,50,1,MPosX1-10,MPosX1+10); MPX1=Pars(1);SigmaPosX1=Pars(2);
%drawnow
%MPosY1=mean(PY1);[Pars]=hfitg(PY1,50,1,MPosY1-10,MPosY1+10); MPY1=Pars(1);SigmaPosY1=Pars(2);
%drawnow

%disp('		Scintilator cuts')
%Aap=(Ch13+Ch14)/2; % comment to run after CAT
MAap1=mean(Aap1);[Pars]=hfitg(Aap1,50,1,MAap1-100,MAap1); MAap1=Pars(1);SigmaAap1=Pars(2);
drawnow
MAap2=mean(Aap2);[Pars]=hfitg(Aap2,50,1,MAap2-100,MAap2); MAap2=Pars(1);SigmaAap2=Pars(2);
drawnow

%disp('		Tracking hit strips cut')
%	abs(HitsX1(:,1)-3)<3&...
%	abs(HitsY1(:,1)-3)<3&...
%	abs(HitsX2(:,1)-3)<3&...
%	abs(HitsY2(:,1)-3)<3&...
%	abs(HitsX3(:,1)-3)<3&...
%	abs(HitsY3(:,1)-3)<3&...

As_=(sort(Arpc));
I=1:110000;
TsMcut=1*50;

[Tc1,Tc2,TsM_,A_,I,Tstart1_]=cut(find(...
	Arpc>As_(length(Arpc)*.02)&...
	Arpc<1000 &...
	abs(Aap1-MAap1)<SigmaAap1*2&...
	abs(Aap2-MAap2)<SigmaAap2*2&...
	abs(TsM)<TsMcut &... 
	abs((HitsX1)-3)<3&...
	abs((HitsY1)-3)<3&...
	abs((HitsX2)-3)<3&...
	abs((HitsY2)-3)<3&...
	abs(x)<5&...
	abs(y)<5&...
	abs(TCorrT1)<1000&...
	abs(TCorrT2)<1000&...
	abs(TsM)<1000&...
	abs(Tstart1-635)<15&...
	1),TCorrT1,TCorrT2,TsM,Arpc,I,Tstart1);

%	abs(PX1-MPosX1)<SigmaPosX1*2&... 
%	abs(PY1-MPosY1)<SigmaPosY1*2&... 
%	TCorr3<5&...

%	abs(HitsX2-3)<3&...
%	abs(HitsY2-3)<3&...

%	Tevent>100 & Tevent<200 &...
%	abs(PosX1-MPosX1)<SigmaPosX1*1&... 
%	abs(PosY1-MPosY1)<SigmaPosY1*1&... 

Events=length(Tc1)

Tc1=[];
Tc2=[];
TsM_=[];

disp('--------realign the distributions------------');

for i=1:10000:100000
disp('-------')
U=I(find(I>i & I<i+10000));
%J=find(abs(TCorrT1(U))<100);%mean(TCorrT1(U(J)))
T1=TCorrT1(U);
T2=TCorrT2(U);
Ts=TsM(U);
A=Arpc(U);
T1=uncorr(T1,Ts);T1=uncorr(T1,A);
[Mean, Sigma, Max]=stdfit(T1);Mean,drawnow;Tc1=[Tc1; T1-Mean];
T2=uncorr(T2,Ts);T2=uncorr(T2,A);
[Mean, Sigma, Max]=stdfit(T2);Mean,drawnow;Tc2=[Tc2; T2-Mean];
%[Mean, Sigma, Max]=stdfit(Ts,num2str(i));Mean,drawnow;
TsM_=[TsM_; Ts];

end

TsM=TsM_;

% calculate sigmas

disp(' ------- RPC, T1,T2 ------------------------------------');
[Mean, Sigma, Max]=stdfit(Tc1,'10');Sigma10=Sigma
[Mean, Sigma, Max]=stdfit(Tc2,'20');Sigma20=Sigma
[N,X]=histf(TsM_,-TsMcut:50:TsMcut);[Pars]=hfitg(X,N,1); Sigma12=Pars(2),title('TsM');

M=[1 0 1;1 1 0;0 1 1];
Sigmas=sqrt(M\([Sigma10;Sigma20;Sigma12].^2))


end % basic processing

%tails
disp(' ------- RPC, T1,T2 with tails ------------------------------------');

par=hfitg(Tc1,50,0,-500,500);logy;yaxis(1,10000);
Mean=par(1);Sigma=par(2);Max=par(3);
Tc1=Tc1-Mean;
par=hfitg(Tc1,50,1,-3*Sigma,3*Sigma);logy;yaxis(1,10000);
Mean=par(1);Sigma10=par(2),Max=par(3);
title('Rpc-T1')

par=hfitg(Tc2,50,0,-500,500);logy;yaxis(1,10000);
Mean=par(1);Sigma=par(2);Max=par(3);
Tc2=Tc2-Mean;
par=hfitg(Tc2,50,1,-3*Sigma,3*Sigma);logy;yaxis(1,10000);
Mean=par(1);Sigma20=par(2),Max=par(3);
title('Rpc-T2')

%par=hfitg(TsM,50,0,-500,500);logy;yaxis(1,10000);
%Mean=par(1);Sigma=par(2);Max=par(3);
%TsM=TsM-Mean;
%par=hfitg(TsM,50,1,-TsMcut*.9,TsMcut*.9);logy;yaxis(1,10000);Mean=par(1);Sigma12=par(2),Max=par(3);
[N,X]=histf(TsM_,-TsMcut:50:TsMcut);[Pars]=hfitg(X,N,1); Sigma12=Pars(2),title('TsM');

title('T1-T2')

M=[1 0 1;1 1 0;0 1 1];
Sigmas=sqrt(M\([Sigma10;Sigma20;Sigma12].^2))

% tails plots

[N1,X1]=histf(Tc1,-500:20:500);
figure;stairs(X1,N1/max(N1)); hold on; xaxis(-500,500)
plot(X1+10,g(X1,Mean,Sigma),'-b');
yaxis(0.011,1);logy;grid on; hold off
title('Rpc-T1')
J=find(N1<max(N1)*0.011);N1_3s=N1;N1_3s(J)=J*0; % 3 sigma
J=find(N1<max(N1)*0.011);N1_2s=N1/max(N1(J));N1_2s=J*0; % 2 sigma

[N,X]=histf(Tc2,-500:20:500);
figure;stairs(X,N/max(N)); hold on; xaxis(-500,500)
plot(X+10,g(X,Mean,Sigma),'-b');
yaxis(0.011,1);logy;grid on; hold off
title('Rpc-T2')

%show all
bins=-1000:20:1000;
figure;[N,X]=histf([Tc1' bins],bins);hfitg(X,N,1);
title('Rpc-T1'); logy;yaxis(1,2000);

return

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

