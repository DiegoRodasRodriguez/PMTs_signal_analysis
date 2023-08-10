close all;
clearvars;

dirQE    = 'E:\HOME_RareEventsGroup\Diego\QE_Calc\QEdata\';
dirYield = 'E:\HOME_RareEventsGroup\Diego\QE_Calc\NobleLightGenerator\';

lambda_ = 100:1:450;

%%%%%%%%%%%%%%%%%%%%%%% LOAD SCINTILLATION FILES %%%%%%%%%%%%%%%%%%%%%%%%%%

file = 'ArThird05bar'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_ArThird05bar = U(:,1); YieldArThird05bar = U(:,2);

file = 'ArThird1bar'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_ArThird1bar = U(:,1); YieldArThird1bar = U(:,2);

file = 'ArThird10bar'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_ArThird10bar = U(:,1); YieldArThird10bar = U(:,2);

file = 'ArThird18bar'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_ArThird18bar = U(:,1); YieldArThird18bar = U(:,2);

file = 'XeThird1bar'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_XeThird1bar = U(:,1); YieldXeThird1bar = U(:,2);

file = 'XeLightKoehler'; %It would be nice to include these as a generator!
eval(['load ', dirYield, file]);
lambda_XeSecond = U(:,1); YieldXeSecond = U(:,2);    %Koehler was taken at 34bar, but it should not matter

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
      
   YieldArThird05bar_(i)  = interp1(lambda_ArThird05bar,   YieldArThird05bar, lambda_(i), 'linear', 0);   
   YieldArThird1bar_(i)   = interp1(lambda_ArThird1bar,    YieldArThird1bar,  lambda_(i), 'linear', 0);
   YieldArThird10bar_(i)  = interp1(lambda_ArThird10bar,   YieldArThird10bar, lambda_(i), 'linear', 0);
   YieldArThird18bar_(i)  = interp1(lambda_ArThird18bar,   YieldArThird18bar, lambda_(i), 'linear', 0);
   
   YieldXeThird1bar_(i)   = interp1(lambda_XeThird1bar,    YieldXeThird1bar,  lambda_(i), 'linear', 'extrap');   
   YieldXeSecond_(i)      = interp1(lambda_XeSecond,       YieldXeSecond,     lambda_(i), 'linear', 0);
   
   T_Filter250to400_(i)   = interp1(lambda_Filter250to400, T_Filter250to400,  lambda_(i), 'linear', 0);
   T_FilterAbove250_(i)   = interp1(lambda_FilterAbove250, T_FilterAbove250,  lambda_(i), 'linear', 0);
   T_FilterVisible_(i)    = interp1(lambda_FilterVisible,  T_FilterVisible,   lambda_(i), 'linear', 0);
     
end
 
QE_5070_          (QE_5070_<0)          = 0;
QE_R7378_         (QE_R7378_<0)         = 0;
YieldXeThird1bar_ (YieldXeThird1bar_<0) = 0;

YieldArThird05bar_ =  YieldArThird05bar_/sum(YieldArThird05bar_);
YieldArThird1bar_  =  YieldArThird1bar_/sum(YieldArThird1bar_);
YieldArThird10bar_ =  YieldArThird10bar_/sum(YieldArThird10bar_);
YieldArThird18bar_ =  YieldArThird18bar_/sum(YieldArThird18bar_);

YieldXeThird1bar_  =  YieldXeThird1bar_/sum(YieldXeThird1bar_);
YieldXeSecond_     =  YieldXeSecond_/sum(YieldXeSecond_);

figure; subplot(2,1,1);
title('Light Yield (argon)'); hold on;
plot(lambda_, YieldArThird05bar_, 'b-');
plot(lambda_, YieldArThird1bar_,  'r-');
plot(lambda_, YieldArThird10bar_, 'g-'); 
plot(lambda_, YieldArThird18bar_, 'k-'); 
box; title('argon');
xlabel('wavelength [nm]');
legend('P=0.5bar', 'P=1bar', 'P=10bar', 'P=18bar');

subplot(2,1,2);
title('Light Yield (xenon)'); hold on;
plot(lambda_, YieldXeThird1bar_, 'b-');
plot(lambda_, YieldXeSecond_,    'r-');
box; title('xenon');
xlabel('wavelength [nm]');
legend('(3rd cont)P=1bar','(2nd cont)P=34bar');

figure; title('QE'); hold on;
plot(lambdaQE_R7378,     QE_R7378,     'b-');
plot(lambdaQE_5070,      QE_5070,      'r-');
plot(lambda_Filter250to400, T_Filter250to400, 'g-');
plot(lambda_FilterAbove250, T_FilterAbove250, 'k-');
plot(lambda_FilterVisible,  T_FilterVisible,  'c-');
box; legend('QE R7378', 'QE 5070', 'T filter 1', 'T filter 2', 'T filter 3');
xlabel('wavelength [nm]');

QE_ArThird05bar_R7378 = sum(YieldArThird05bar_ .* QE_R7378_)
QE_ArThird1bar_R7378  = sum(YieldArThird1bar_  .* QE_R7378_)
QE_ArThird10bar_R7378 = sum(YieldArThird10bar_ .* QE_R7378_)
QE_ArThird18bar_R7378 = sum(YieldArThird18bar_ .* QE_R7378_)
QE_XeThird1bar_R7378  = sum(YieldXeThird1bar_  .* QE_R7378_)
QE_XeSecond_R7378     = sum(YieldXeSecond_     .* QE_R7378_)

QE_ArThird05bar_5070 = sum(YieldArThird05bar_ .* QE_5070_)
QE_ArThird1bar_5070  = sum(YieldArThird1bar_  .* QE_5070_)
QE_ArThird10bar_5070 = sum(YieldArThird10bar_ .* QE_5070_)
QE_ArThird18bar_5070 = sum(YieldArThird18bar_ .* QE_5070_)
QE_XeThird1bar_5070  = sum(YieldXeThird1bar_  .* QE_5070_)

QE_XeSecond_5070     = sum(YieldXeSecond_     .* QE_5070_)

QE_ArThird05bar_R7378_filter250to400 = sum(YieldArThird05bar_ .* QE_R7378_ .* T_Filter250to400_)
QE_ArThird1bar_R7378_filter250to400  = sum(YieldArThird1bar_  .* QE_R7378_ .* T_Filter250to400_)
QE_ArThird10bar_R7378_filter250to400 = sum(YieldArThird10bar_ .* QE_R7378_ .* T_Filter250to400_)
QE_ArThird18bar_R7378_filter250to400 = sum(YieldArThird18bar_ .* QE_R7378_ .* T_Filter250to400_)
QE_XeThird1bar_R7378_filter250to400  = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_Filter250to400_)

QE_XeThird1bar_R7378_filter250to400 = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_Filter250to400_)

QE_XeThird1bar_R7378_filterAbove250 = sum(YieldXeThird1bar_  .* QE_R7378_ .* T_FilterAbove250_)
QE_XeThird1bar_5070_filterVisible   = sum(YieldXeThird1bar_  .* QE_5070_ .* T_FilterVisible_)


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