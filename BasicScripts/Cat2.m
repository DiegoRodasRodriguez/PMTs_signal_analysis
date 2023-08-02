% concatenate runs (after TA corr)

clear all
close all

TcorrT1_=[];
TcorrT2_=[];
TsM_=[];
Tstart1_=[];
SigmaRpcT1_=[];
SigmaRpcT2_=[];
SigmaT0T1_=[];
Arpc_=[];
Aap1_=[];
Aap2_=[];
Tevent_=[];
Baseline_=[];
HitsX1_=[];
HitsY1_=[];
HitsX2_=[];
HitsY2_=[];
PosX1_=[];
PosY1_=[]; 
PosX2_=[];
PosY2_=[];


for run=[2108 2112:2116 2118:2122]
disp('-------------------');
s=['load data\corr' num2str(run)]
eval(s)

SigmaRpcT1=SigmaRpcT1
SigmaRpcT2=SigmaRpcT2
SigmaT0T1=SigmaT0T1

%I=find(abs(TcorrT1)<300);TcorrT1=TcorrT1-mean(TcorrT1(I));
%I=find(abs(TcorrT2)<300);TcorrT2=TcorrT2-mean(TcorrT2(I));
I=find(abs(TsM)<6);TsM=TsM-mean(TsM(I));

TCorrT1_=[TCorrT1_;TcorrT1];
TCorrT2_=[TCorrT2_;TcorrT2];
TsM_=[TsM_;TsM];
Tstart1_=[Tstart1_;Tstart1];
Arpc_=[Arpc_;Arpc];
Aap1_=[Aap1_;Aap1];
Aap2_=[Aap2_;Aap2];
HitsX1_=[HitsX1_;HitsX1(:,1)];
HitsY1_=[HitsY1_;HitsY1(:,1)];
HitsX2_=[HitsX2_;HitsX2(:,1)];
HitsY2_=[HitsY2_;HitsY2(:,1)];
PX1=[PX1;PosX1(:,1)];
PY1=[PY1;PosY1(:,1)];
PX2=[PX2;PosX2(:,1)];
PY2=[PY2;PosY2(:,1)];

end

TCorrT1=TCorrT1_;
TCorrT2=TCorrT2_;
TsM=TsM_;
Tstart1=Tstart1_;
Arpc=Arpc_;
Aap1=Aap1_;
Aap2=Aap2_;
HitsX1=HitsX1_;
HitsY1=HitsY1_;
HitsX2=HitsX2_;
HitsY2=HitsY2_;

s=['save data\CorrAll2 ' 'TCorrT1 TCorrT2 TsM Tstart1 Arpc Aap1 Aap2 HitsX1 HitsY1 HitsX2 HitsY2 PX1 PY1 PX2 PY2 '];
eval(s)
