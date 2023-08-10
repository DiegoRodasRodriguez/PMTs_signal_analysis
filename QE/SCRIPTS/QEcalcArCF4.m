close all;
clearvars;

dirQE    = 'E:\HOME_RareEventsGroup\Diego\QE_Calc\QEdata\';
dirYield = 'E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\';

lambda_ = 100:1:900;

Yield90Ar_ = LightSpectrumInArCF4_sec(0.1, lambda_);
Yield95Ar_ = LightSpectrumInArCF4_sec(0.05, lambda_);

Yield90Ar_UV   = Yield90Ar_(find(lambda_>250 & lambda_<400));
lambda90Ar_UV  = lambda_   (find(lambda_>250 & lambda_<400));
Yield90Ar_vis  = Yield90Ar_(find(lambda_>400 & lambda_<740));
lambda90Ar_vis = lambda_   (find(lambda_>400 & lambda_<740));
Yield95Ar_UV   = Yield95Ar_(find(lambda_>250 & lambda_<400));
lambda95Ar_UV  = lambda_   (find(lambda_>250 & lambda_<400));
Yield95Ar_vis  = Yield95Ar_(find(lambda_>400 & lambda_<740));
lambda95Ar_vis = lambda_   (find(lambda_>400 & lambda_<740));

%%%%%%%%%%%%%%%%%%%%% LOAD QE/TRANSMISSION FILES %%%%%%%%%%%%%%%%%%%%%%%%%%

file = 'QE_5070';
eval(['load ', dirQE, file]);
lambdaQE_5070  = U(:,1); QE_5070   = U(:,2)/100;

file = 'R7378_QE';
eval(['load ', dirQE, file]);
lambdaQE_R7378 = U(:,1); QE_R7378  = U(:,2);

file = 'Filter250to400';
eval(['load ', dirQE, file]);
lambda_Filter250to400 = U(:,1); T_Filter250to400  = U(:,2)/100;

file = 'FilterAbove250';
eval(['load ', dirQE, file]);
lambda_FilterAbove250 = U(:,1); T_FilterAbove250  = U(:,2)/100;

file = 'FilterVisible';
eval(['load ', dirQE, file]);
lambda_FilterVisible = U(:,1); T_FilterVisible  = U(:,2)/100;


for i=1:length(lambda_) 
      
   QE_5070_(i)            = interp1(lambdaQE_5070,         QE_5070,           lambda_(i), 'linear', 'extrap');    
   QE_R7378_(i)           = interp1(lambdaQE_R7378,        QE_R7378,          lambda_(i), 'linear', 'extrap'); 
       
   T_Filter250to400_(i)   = interp1(lambda_Filter250to400, T_Filter250to400,  lambda_(i), 'linear', 0);
   T_FilterAbove250_(i)   = interp1(lambda_FilterAbove250, T_FilterAbove250,  lambda_(i), 'linear', 0);
   T_FilterVisible_(i)    = interp1(lambda_FilterVisible,  T_FilterVisible,   lambda_(i), 'linear', 0);
      
   Yield90Ar_UV_(i)   = interp1(lambda90Ar_UV,  Yield90Ar_UV,  lambda_(i), 'linear', 0);
   Yield90Ar_vis_(i)  = interp1(lambda90Ar_vis, Yield90Ar_vis, lambda_(i), 'linear', 0);
   Yield95Ar_UV_(i)   = interp1(lambda95Ar_UV,  Yield95Ar_UV,  lambda_(i), 'linear', 0);
   Yield95Ar_vis_(i)  = interp1(lambda95Ar_vis, Yield95Ar_vis, lambda_(i), 'linear', 0);
   
end
 
QE_5070_          (QE_5070_<0)          = 0;
QE_R7378_         (QE_R7378_<0)         = 0;

Yield90Ar_UV_  = Yield90Ar_UV_/sum(Yield90Ar_UV_);
Yield95Ar_UV_  = Yield95Ar_UV_/sum(Yield95Ar_UV_);
Yield90Ar_vis_ = Yield90Ar_vis_/sum(Yield90Ar_vis_);
Yield95Ar_vis_ = Yield95Ar_vis_/sum(Yield95Ar_vis_);

figure; subplot(2,1,1);
title('Light Yield (argon + CF4)'); hold on;
plot(lambda_, Yield90Ar_,  'b-');
plot(lambda_, Yield95Ar_,  'r-');
box; title('argon');
xlabel('wavelength [nm]');
legend('f_{CF4}=10%', 'f_{CF4}=5%');

figure; title('QE'); hold on;
plot(lambdaQE_R7378,     QE_R7378,     'b-');
plot(lambdaQE_5070,      QE_5070,      'r-');
plot(lambda_Filter250to400, T_Filter250to400, 'g-');
plot(lambda_FilterAbove250, T_FilterAbove250, 'k-');
plot(lambda_FilterVisible,  T_FilterVisible,  'c-');
box; legend('QE R7378', 'QE 5070', 'T filter 1', 'T filter 2', 'T filter 3');
xlabel('wavelength [nm]');

% 
QE_90Ar_UV_R7378  = sum(Yield90Ar_UV_  .* QE_R7378_ .* T_Filter250to400_)
QE_90Ar_vis_R7378 = sum(Yield90Ar_vis_ .* QE_R7378_ .* T_FilterVisible_)
QE_95Ar_UV_R7378  = sum(Yield95Ar_UV_  .* QE_R7378_ .* T_Filter250to400_)
QE_95Ar_vis_R7378 = sum(Yield95Ar_vis_  .* QE_R7378_ .* T_FilterVisible_)

QE_90Ar_UV_5070   = sum(Yield90Ar_UV_  .* QE_5070_  .* T_Filter250to400_)
QE_90Ar_vis_5070  = sum(Yield90Ar_vis_ .* QE_5070_  .* T_FilterVisible_)
QE_95Ar_UV_5070   = sum(Yield95Ar_UV_  .* QE_5070_  .* T_Filter250to400_)
QE_95Ar_vis_5070  = sum(Yield95Ar_vis_  .* QE_5070_  .* T_FilterVisible_)

% QE_ArThird1bar_R7378  = sum(YieldArThird1bar_  .* QE_R7378_)
% QE_ArThird10bar_R7378 = sum(YieldArThird10bar_ .* QE_R7378_)
% QE_ArThird18bar_R7378 = sum(YieldArThird18bar_ .* QE_R7378_)
% QE_XeThird1bar_R7378  = sum(YieldXeThird1bar_  .* QE_R7378_)
% QE_XeSecond_R7378     = sum(YieldXeSecond_     .* QE_R7378_)
% 
% QE_ArThird05bar_5070 = sum(YieldArThird05bar_ .* QE_5070_)
% QE_ArThird1bar_5070  = sum(YieldArThird1bar_  .* QE_5070_)
% QE_ArThird10bar_5070 = sum(YieldArThird10bar_ .* QE_5070_)
% QE_ArThird18bar_5070 = sum(YieldArThird18bar_ .* QE_5070_)
% QE_XeThird1bar_5070  = sum(YieldXeThird1bar_  .* QE_5070_)
% 
% QE_XeSecond_5070     = sum(YieldXeSecond_     .* QE_5070_)
% 
% QE_ArThird05bar_R7378_filter250to400 = sum(YieldArThird05bar_ .* QE_R7378_ .* T_Filter250to400_)
% QE_ArThird1bar_R7378_filter250to400  = sum(YieldArThird1bar_  .* QE_R7378_ .* T_Filter250to400_)
% QE_ArThird10bar_R7378_filter250to400 = sum(YieldArThird10bar_ .* QE_R7378_ .* T_Filter250to400_)
% QE_ArThird18bar_R7378_filter250to400 = sum(YieldArThird18bar_ .* QE_R7378_ .* T_Filter250to400_)
% QE_XeThird1bar_R7378_filter250to400  = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_Filter250to400_)
% 
% QE_XeThird1bar_R7378_filter250to400 = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_Filter250to400_)
% 
% QE_XeThird1bar_R7378_filterAbove250 = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_FilterAbove250_)
% QE_XeThird1bar_5070_filterVisible   = sum(YieldXeThird1bar_  .* QE_5070_ .* T_FilterVisible_)


% QE_95Ar_R1926   = sum(Yield95Ar_   .* Transm_ .* QE_R1926_)
% QE_60He_R1926   = sum(Yield60He_   .* Transm_ .* QE_R1926_)
% QE_100CF4_R1926 = sum(Yield100CF4_ .* Transm_ .* QE_R1926_)
% 
% QE_67Ar_R7378   = sum(Yield67Ar_   .* Transm_ .* QE_R7378_)
% QE_95Ar_R7378   = sum(Yield95Ar_   .* Transm_ .* QE_R7378_)
% QE_60He_R7378   = sum(Yield60He_   .* Transm_ .* QE_R7378_)
% QE_100CF4_R7378 = sum(Yield100CF4_ .* Transm_ .* QE_R7378_)
% 
% QE_67Ar_5070    = sum(Yield67Ar_   .* Transm_ .* QE_5070_)
% QE_95Ar_5070    = sum(Yield95Ar_   .* Transm_ .* QE_5070_)
% QE_60He_5070    = sum(Yield60He_   .* Transm_ .* QE_5070_)
% QE_100CF4_5070  = sum(Yield100CF4_ .* Transm_ .* QE_5070_)
% 


mosaic;