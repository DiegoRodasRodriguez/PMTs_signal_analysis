clear all;
close all;
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_ArCF4_90_10_1bar_GEM_UVcut;
LambdaMin = 400; LambdaMax = 870; LambdaStep = 1;
LambdaMin_ = LambdaMin - 50;
LambdaMax_ = LambdaMax + 50;

Lambda_ = LambdaMin_:LambdaStep:LambdaMax_;

Lambda = U(:,1);
Yield  = U(:,2)*1e-3; %ph/el/nm

Yield_ = interp1(Lambda, Yield, Lambda_, 'linear', 'extrap');    
Yield_(Yield_<0)=0;
Norm   = sum(Yield_(Lambda_>400));
Yield_ = Yield_/Norm * 0.55;
Yield_(Lambda_<LambdaMin | Lambda_>LambdaMax) = -1e-6;

figure;
plot(Lambda, Yield, 'r');
figure;
plot(Lambda_, Yield_, '-');

Yield  = Yield_;
Lambda = Lambda_;

%save E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_ArCF4_90_10_1bar_GEM_UVcut_NORM ...
%    Lambda Yield;