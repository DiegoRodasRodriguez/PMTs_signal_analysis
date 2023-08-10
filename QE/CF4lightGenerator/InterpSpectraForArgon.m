close all;
clear all;

LambdaMax  = 900;
LambdaMin  = 200;
LambdaStep = 1;
Lambda_    = LambdaMin:LambdaStep:LambdaMax;
NobleGas      = 'Ar';
AdmixtureCF4  = 0.2;
Admixture_    = 0:0.05:1;
YieldInterp   = zeros(length(Admixture_), length(Lambda_));
phEl_UV       = zeros(size(Admixture_));
phEl_vis      = zeros(size(Admixture_));
phEl_IR       = zeros(size(Admixture_));

load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_PureCF4_1bar_GEM_NORM;
LambdaPure = Lambda; YieldPure = Yield;
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_HeCF4_60_40_1bar_GEM_NORM;
LambdaHe60 = Lambda; YieldHe60 = Yield; 
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_ArCF4_95_5_1bar_GEM_NORM;
LambdaAr95 = Lambda; YieldAr95 = Yield;
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_ArCF4_90_10_1bar_GEM_UVcut_NORM;
LambdaAr90 = Lambda; YieldAr90 = Yield;
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_ArCF4_33_67_1bar_GEM_NORM;
LambdaAr33 = Lambda; YieldAr33 = Yield;
load E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\CF4lightGenerator\S2_PureAr_1bar_GEM_NORM;
LambdaPureAr = Lambda; YieldPureAr = Yield;

figure; hold on;
plot(LambdaPure, YieldPure, '-k');
plot(LambdaAr33, YieldAr33, '-r');
plot(LambdaAr90, YieldAr90, '-g');
plot(LambdaAr95, YieldAr95, '-b');

legend('pure CF_4','Ar/CF_4, 33/67', 'Ar/CF_4, 90/10', 'Ar/CF_4, 95/5');

figure; hold on;
plot(LambdaPure, YieldPure, '-k');
plot(LambdaHe60, YieldHe60, '-r');

sum(YieldAr33(LambdaAr33>400))
sum(YieldAr90(LambdaAr90>400))
sum(YieldAr95(LambdaAr95>400))
sum(YieldHe60(LambdaHe60>400))
sum(YieldPure(LambdaPure>400))

for i=1:length(Admixture_), YieldInterp(i,:) = LightSpectrumInArCF4(Admixture_(i), Lambda_); end

figure;
subplot(2,5,1);  plot(Lambda_, YieldInterp(1,:),  '-');
subplot(2,5,2);  plot(Lambda_, YieldInterp(2,:),  '-');
subplot(2,5,3);  plot(Lambda_, YieldInterp(3,:),  '-');
subplot(2,5,4);  plot(Lambda_, YieldInterp(7,:),  '-');
subplot(2,5,5);  plot(Lambda_, YieldInterp(9,:),  '-');
subplot(2,5,6);  plot(Lambda_, YieldInterp(11,:),  '-');
subplot(2,5,7);  plot(Lambda_, YieldInterp(13,:),  '-');
subplot(2,5,8);  plot(Lambda_, YieldInterp(15,:),  '-');
subplot(2,5,9);  plot(Lambda_, YieldInterp(18,:),  '-');
subplot(2,5,10); plot(Lambda_, YieldInterp(21,:),  '-');

for i=1:length(Admixture_)
    YieldTmp    = YieldInterp(i,:);
    phEl_UV(i)  = sum(YieldTmp(Lambda_>=LambdaMin & Lambda_<=400));
    phEl_vis(i) = sum(YieldTmp(Lambda_>=400       & Lambda_<=700));
    phEl_IR(i)  = sum(YieldTmp(Lambda_>=700       & Lambda_<=LambdaMax));
end

%phEl_TOT = phEl

%Morozov and Fraga together with errors
%Morozov: 60/50/140
%Fraga: 80/60/140 o 60/40/140 (unclear) NOTE: indeed it is not clear which one.
phEl_TOT_LIP     = [0.575, 0.55, 0.3, 0.09];
Err_phEl_TOT_LIP = [0.575*0.2, 0.55*0.2, 0.3*0.2, 0.09*sqrt(0.17^2+0.10^2)];
conc_TOT_LIP     = [0.05, 0.1, 0.67, 1];

%Fraga with CCD camera (2GEMs, hole/kapton/pitch = 80/50/140
NORM                 = (0.6*1e-5)/0.09;
phEl_TOT_LIP_CCD     = [2.2, 1.9, 1.7, 1.2, 0.8]*1e-5/NORM;
Err_phEl_TOT_LIP_CCD = phEl_TOT_LIP_CCD*0.2; %FIXME: error CCD?
conc_TOT_LIP_CCD     = [0.05, 0.2, 0.4, 0.8, 1];

%Note: 80/60/140 can give up to a factor x2 less.


figure; hold on;
plot(Admixture_(Admixture_>=0.05), phEl_UV(Admixture_>=0.05), '-');
plot(Admixture_(Admixture_>=0.05), phEl_vis(Admixture_>=0.05),'-r');
plot(Admixture_(Admixture_>=0.05), phEl_IR(Admixture_>=0.05), '-g');
plot(Admixture_(Admixture_>=0.05), phEl_vis(Admixture_>=0.05)+phEl_IR(Admixture_>=0.05), '-k');
errorbar(conc_TOT_LIP, phEl_TOT_LIP, Err_phEl_TOT_LIP, 's');
errorbar(conc_TOT_LIP_CCD, phEl_TOT_LIP_CCD, Err_phEl_TOT_LIP_CCD, 'o');
legend('UV', 'vis', 'IR', 'vis+IR', 'single GEM (60um holes)', 'double GEM (80um holes)');

plot(Admixture_(Admixture_<=0.05), phEl_UV(Admixture_<=0.05), ':b');
plot(Admixture_(Admixture_<=0.05), phEl_vis(Admixture_<=0.05),':r');
plot(Admixture_(Admixture_<=0.05), phEl_IR(Admixture_<=0.05), ':g');
plot(Admixture_(Admixture_<=0.05), phEl_vis(Admixture_<=0.05)+phEl_IR(Admixture_<=0.05), ':k');

