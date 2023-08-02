% August 2000 runs (with 2 start counters)

disp('====================================')
disp('== Analyse the sharing of charge  ==')
filename=filename
disp('====================================')

save tmp filename
close all		% comment to use StartStp or mana1
clear all

OFFSETrpc=1;SIGMArpc=1;

TchLeft=17; TchRight=18;AchLeft=5; AchRight=6;

load data\dat_0204;
ped1=(U(:,AchLeft)+U(:,AchRight))/2;
ped2=(U(:,7)+U(:,8))/2;

load tmp
s=['load data\' filename], eval(s); load tmp
%load data\dat_0134;
%load data\350_351;filename='350_351';

if 0
   disp('		Number of events<=20000')
   U=U(1:min(20000,length(U)),:);
end

U(:,[13,14,15,16,TchLeft,TchRight])=U(:,[13,14,15,16,TchLeft,TchRight])+rand(length(U),6); % shake data

Trpc=(U(:,TchLeft)+U(:,TchRight))/2; Arpc=(U(:,AchLeft)+U(:,AchRight))/2;
DTrpc=U(:,TchLeft)-U(:,TchRight); 
Arpc2=(U(:,7)+U(:,8))/2;

Tstart1=(U(:,13)+U(:,14))/2; Tstart2=(U(:,15)+U(:,16))/2; 
DT1=(U(:,13)-U(:,14)); DT2=(U(:,15)-U(:,16));
Aap1=(U(:,1)+U(:,2))/2; Aap2=(U(:,3)+U(:,4))/2; 
Events=length(Arpc);

TrpcT1=Trpc-Tstart1+mean(Tstart1);
TrpcT2=Trpc-Tstart2+mean(Tstart2);
%TrpcT1=Trpc-(Tstart1+Tstart2)/2+mean((Tstart1+Tstart2)/2); for single start counter analysis

% center some distributions

TsM=Tstart1-Tstart2;
[N,X]=histf(TsM,-100:400); TsM=TsM-max(X(find(N==max(N)))); % center peak on 0
[N,X]=histf(TsM,-10:10); TsM=TsM-max(X(find(N==max(N))));

DTrpc=DTrpc-mean(DTrpc); [N,X]=histf(DTrpc,-500:10:500); DTrpc=DTrpc-max(X(find(N==max(N))));% center peak on 0
DT1=DT1-mean(DT1); [N,X]=histf(DT1,-30:30); DT1=DT1-max(X(find(N==max(N))));% center peak on 0
DT2=DT2-mean(DT2); [N,X]=histf(DT2,-30:30); DT2=DT2-max(X(find(N==max(N))));% center peak on 0

%TrpcT1=DT-mean(DT)+500; %AMPLITUDE-CORRECT DT

%disp(' ----------- Efficiency---------------------');

Iok=1:Events;	% master index
Events=length(Iok),IinEff=find(U(Iok,TchLeft)>=(2048) & U(Iok,TchRight)>=2048);
Eff_TOF_Hardware=(1-length(IinEff)/Events)*100
%BadTriggers=length(find(Tstart1>=2048 | Tstart2>=2048))/Events*100


disp(' ------- apply EXTERNAL cuts ---------------------');

OKflag=ones(Events,1); Icut=[]; cutaux;

if 0
disp('		valid Tstart')
	Icut=find((Tstart1>=2048)&(Tstart2>=2048));
	cutaux
end

if 0
   t='		TsM cut';   disp(t)
	Icut=find(abs(TsM)>3);
   cutaux
   figure;histf(TsM,-20:20);
	hold on;histf(TsM(LocalIok),-20:20);lastline('color','b')
	title(t)   
end

if 0
   t='		Scintilator position (DT1,2) cut';   disp(t)
   
   CutSigma=2
   [N,X]=histf(DT1,-10:10); DT1=DT1-max(X(find(N==max(N))));
   [N,X]=histf(DT1,-10:10); 
   Llim=X(max(find(N/max(N)<g(CutSigma,0,1) & X<0)));
   Ulim=X(min(find(N/max(N)<g(CutSigma,0,1) & X>0)));
   Icut=find(DT1<Llim | DT1>Ulim);
   cutaux
   figure;
   subplot(1,2,1);histf(DT1,-20:20);
   hold on;histf(DT1(LocalIok),-20:20);lastline('color','b')
   
   [N,X]=histf(DT2,-10:10); DT2=DT2-max(X(find(N==max(N))));
   Llim=X(max(find(N/max(N)<g(CutSigma,0,1)  & X<0)));
   Ulim=X(min(find(N/max(N)<g(CutSigma,0,1)  & X>0)));
   Icut=find(DT2<Llim | DT2>Ulim);
   cutaux
   subplot(1,2,2);histf(DT2,-20:20);
	hold on;histf(DT2(LocalIok),-20:20);lastline('color','b')
	title(t)   

end

if 0
   t='	Scintilator Amplitude cuts'; disp(t)
   
   CutSigma=2
   [N,X]=histf(Aap1,0:10:2000); Aap1=Aap1-max(X(find(N==max(N))));
   [N,X]=histf(Aap1,-400:10:400); Aap1=Aap1-max(X(find(N==max(N))));
   [N,X]=histf(Aap1,-400:10:400); 
   Llim=X(max(find(N/max(N)<g(CutSigma,0,1) & X<0)));
   Ulim=X(min(find(N/max(N)<g(CutSigma,0,1) & X>0)));
   
   Icut=find(Aap1<Llim | Aap1>Ulim);
   cutaux
   figure;
   subplot(2,1,1); histf(Aap1,-400:10:400);
	hold on;histf(Aap1(LocalIok),-400:10:400);lastline('color','b')
   
   [N,X]=histf(Aap2,0:10:2000); Aap2=Aap2-max(X(find(N==max(N))));
   [N,X]=histf(Aap2,-400:10:400); Aap2=Aap2-max(X(find(N==max(N))));
   Llim=X(max(find(N/max(N)<g(CutSigma,0,1) & X<0)));
   Ulim=X(min(find(N/max(N)<g(CutSigma,0,1) & X>0)));
   
   Icut=find(Aap2<Llim | Aap2>Ulim);

   cutaux
   subplot(2,1,2);histf(Aap2,-400:10:400);
	hold on;histf(Aap2(LocalIok),-400:10:400);lastline('color','b')
	title(t)   
end


TimeEff1=Eff

IinEff=find(OKflag & U(:,19)>=(2048) & U(:,20)>=2048);
TimeEff2=(1-length(IinEff)/Events)*100

TrpcT1=TrpcT1(Iok);
TrpcT2=TrpcT2(Iok);
DTrpc=DTrpc(Iok);
Arpc=Arpc(Iok);
Arpc2=Arpc2(Iok);
TsM=TsM(Iok);
Aap1=Aap1(Iok);
Aap2=Aap2(Iok);
DT1=DT1(Iok);
DT2=DT2(Iok);
Tstart1=Tstart1(Iok);
Tstart2=Tstart2(Iok);
Trpc=Trpc(Iok);

	t='-------------------- Charge sharing';disp(t)
%   CutSigma=2
   
   % Strip1
	Llim=min(ped1); Ulim=max(ped1); 

%   Icut=find(Arpc<Llim | Arpc>Ulim);
%   LocalOKflag=ones(length(Arpc),1);	LocalOKflag(Icut)=Icut*0;	LocalIok=find(LocalOKflag);

   I1=find(Arpc>Ulim);
   figure; subplot(2,1,1)
   histf(Arpc,0:300); title(t)
%   hold on;histf(Arpc(LocalIok),0:300);lastline('color','b')
   hold on;histf(Arpc(I1),0:300);lastline('color','b')
   histf(ped1,0:300),lastline('color','y')
   title(t)   
   N1Good=length(I1)/Events
   
   % Strip2
	Llim2=min(ped2); Ulim2=max(ped2); 
   
%   Icut=find(Arpc2<Llim2 | Arpc2>Ulim2);
%	LocalOKflag=ones(length(Arpc),1);	LocalOKflag(Icut)=Icut*0;	LocalIok=find(LocalOKflag);
   
   I2=find(Arpc2>Ulim2);
   subplot(2,1,2);
   histf(Arpc2,0:300); title(t)
   hold on;histf(Arpc2(I2),0:300);lastline('color','b')
   histf(ped2,0:300),lastline('color','y')
   N2Good=length(I2)/Events

	Igood=find(Arpc>Ulim & Arpc2>Ulim2);
	N12Good=length(Igood)/Events
   N12NoGood=length(find(Arpc<Ulim & Arpc2<Ulim2))/Events
   
   if strcmp(filename,'dat_0147') & Events>20000
      IGood1=Igood(find(Arpc(Igood)>Arpc2(Igood)));
      IGood2=Igood(find(Arpc(Igood)<=Arpc2(Igood)));
   end

   return
   

