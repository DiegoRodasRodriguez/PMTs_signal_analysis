close all;
clear all;

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