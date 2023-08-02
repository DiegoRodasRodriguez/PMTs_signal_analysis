% concatenate runs (after TA corr)

clear all
close all

Tevent_=[];
Tstart_=[];
ARpc1_=[];
ARpc2_=[];
ARpc3_=[];
TCorr1_=[];
TCorr2_=[];
TCorr3_=[];
Aap_=[];
HitsX1_=[];
HitsY1_=[];
HitsX2_=[];
HitsY2_=[];
PosX1_=[];
PosY1_=[];

for run=1434:1439
s=['load data\corr' num2str(run)]
eval(s)

Tevent_=[Tevent_;Tevent];
%Tstart_=[Tstart_;Tstart];
ARpc1_=[ARpc1_;ARpc1];
ARpc2_=[ARpc2_;ARpc2];
ARpc3_=[ARpc3_;ARpc3];
TCorr1_=[TCorr1_;TCorr1];
TCorr2_=[TCorr2_;TCorr2];
TCorr3_=[TCorr3_;TCorr3];
Aap_=[Aap_;Aap];
HitsX1_=[HitsX1_;HitsX1(:,1)];
HitsY1_=[HitsY1_;HitsY1(:,1)];
HitsX2_=[HitsX2_;HitsX2(:,1)];
HitsY2_=[HitsY2_;HitsY2(:,1)];
PosX1_=[PosX1_;PosX1(:,1)];
PosY1_=[PosY1_;PosY1(:,1)];

end

Tevent=Tevent_;
%Tstart=Tstart_;
ARpc1=ARpc1_;
ARpc2=ARpc2_;
ARpc3=ARpc3_;
TCorr1=TCorr1_;
TCorr2=TCorr2_;
TCorr3=TCorr3_;
Aap=Aap_;
HitsX1=HitsX1_;
HitsY1=HitsY1_;
HitsX2=HitsX2_;
HitsY2=HitsY2_;
PosX1=PosX1_;
PosY1=PosY1_;


s=['save data\CorrAll1' ' TCorr1 TCorr2 TCorr3 ' ...
'ARpc1 ARpc2 ARpc3 Aap Tevent HitsX1 HitsY1 HitsX2 HitsY2 PosX1 PosY1'];
eval(s)
