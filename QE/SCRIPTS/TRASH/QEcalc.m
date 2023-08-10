close all;
clearvars;

lambdaCut   = 200;
npr         = 300;
dGem2Plexi  = 0.2;%*0.3/0.2;    %<0.3
tPlexi      = 0.3;%*0.5/0.3;
dPlexiMesh  = 0.25;%*0.5/0.25;
dMeshPMT    = 0.25;%*0.5/0.25;
dPlexiUp2PMT = 1.9;

GainGEM     = 10000;
GainLastGEM = 2*(GainGEM)^0.3; %approx
npe         = 80;

dir  = 'E:\HOME_RareEventsGroup\SETUPS\BAT\ANALYSIS\QE_Calc\CERN_QEdata\';

lambda_ = 200:1:900;

Yield67Ar_  = LightSpectrumInArCF4(0.33, lambda_);
Yield60Ar_  = LightSpectrumInArCF4(0.4, lambda_);
Yield70Ar_  = LightSpectrumInArCF4(0.3, lambda_);
Yield80Ar_  = LightSpectrumInArCF4(0.2, lambda_);
Yield90Ar_  = LightSpectrumInArCF4(0.1, lambda_);
Yield95Ar_  = LightSpectrumInArCF4(0.05, lambda_);
Yield98Ar_  = LightSpectrumInArCF4(0.02, lambda_);
Yield60He_  = LightSpectrumInHeCF4(0.4, lambda_);

%Yield that might be good to include!

file = 'CF4_scintillation';
eval(['load ', dir, file]);
lambda_CF4_tmp = U(:,1); Yield100CF4 = U(:,2);

Yield60He_  = LightSpectrumInHeCF4(0.4, lambda_);

%NOT TOO BAD PMTs

% file = 'QE_H7826';
% eval(['load ', dir, file]);
% lambdaQE_H7826 = U(:,1); QE_H7826  = U(:,2)/100;

file = 'QE_H7422_40';
eval(['load ', dir, file]);
lambdaQE_H7422_40 = U(:,1); QE_H7422_40  = U(:,2)/100;

file = 'QE_H10721_01';
eval(['load ', dir, file]);
lambdaQE_H10721_01 = U(:,1); QE_H10721_01  = U(:,2)/100;

file = 'QE_H10721_20';
eval(['load ', dir, file]);
lambdaQE_H10721_20 = U(:,1); QE_H10721_20  = U(:,2)/100;

file = 'QE_H11461_01';
eval(['load ', dir, file]);
lambdaQE_H11461_01 = U(:,1); QE_H11461_01  = U(:,2)/100;

file = 'QE_H11461_03';
eval(['load ', dir, file]);
lambdaQE_H11461_03 = U(:,1); QE_H11461_03  = U(:,2)/100;

file = 'QE_R375';
eval(['load ', dir, file]);
lambdaQE_R375 = U(:,1); QE_R375  = U(:,2)/100;

%MOST PROMISING PMTs

file = 'QE_PMT_NEXT';
eval(['load ', dir, file]);
lambdaQE_PMT_NEXT = U(:,1); QE_PMT_NEXT  = U(:,2)/100;

file = 'QE_CCD_OrcaFlash';
eval(['load ', dir, file]);
lambdaQE_CCD = U(:,1); QE_CCD  = U(:,2)/100;

file = 'S_H10426';
eval(['load ', dir, file]);
lambdaQE_H10426 = U(:,1); QE_H10426  = U(:,2)/1000*1240./lambdaQE_H10426;

file = 'S_H10426_01';
eval(['load ', dir, file]);
lambdaQE_H10426_01 = U(:,1); QE_H10426_01  = U(:,2)/1000*1240./lambdaQE_H10426_01;

file = 'S_H10425';
eval(['load ', dir, file]);
lambdaQE_H10425 = U(:,1); QE_H10425  = U(:,2)/1000*1240./lambdaQE_H10425;

file = 'S_H10425_01';
eval(['load ', dir, file]);
lambdaQE_H10425_01 = U(:,1); QE_H10425_01  = U(:,2)/1000*1240./lambdaQE_H10425_01;

file = 'S_H7826';
eval(['load ', dir, file]);
lambdaQE_H7826 = U(:,1); QE_H7826  = U(:,2)/1000*1240./lambdaQE_H7826;

file = 'S_H7826_01';
eval(['load ', dir, file]);
lambdaQE_H7826_01 = U(:,1); QE_H7826_01  = U(:,2)/1000*1240./lambdaQE_H7826_01;

%Nausicaa PMTs
file = 'R1926_QE';
eval(['load ', dir, file]);
lambdaQE_R1926 = U(:,1); QE_R1926  = U(:,2);

file = 'R7378_QE';
eval(['load ', dir, file]);
lambdaQE_R7378 = U(:,1); QE_R7378  = U(:,2);

%TransmissionPlexiGlass
% file = 'TransmPlexiGlass';
% eval(['load ', dir, file]);
% lambdaQE_Transm = U(:,1); Transm  = U(:,2)/100 -0.05;

%TransmissionBorosilicate
file = 'TransmBorosilicate';
eval(['load ', dir, file]);
lambdaQE_Transm = U(:,1)*1000; Transm  = U(:,2) -0.05;

for i=1:length(lambda_)
   QE_H7422_40_(i)  = interp1(lambdaQE_H7422_40,  QE_H7422_40,  lambda_(i), 'linear', 'extrap');
   QE_H10721_01_(i) = interp1(lambdaQE_H10721_01, QE_H10721_01, lambda_(i), 'linear', 'extrap');
   QE_H10721_20_(i) = interp1(lambdaQE_H10721_20, QE_H10721_20, lambda_(i), 'linear', 'extrap');
   QE_H11461_01_(i) = interp1(lambdaQE_H11461_01, QE_H11461_01, lambda_(i), 'linear', 'extrap');
   QE_H11461_03_(i) = interp1(lambdaQE_H11461_03, QE_H11461_03, lambda_(i), 'linear', 'extrap');
   QE_R375_(i)      = interp1(lambdaQE_R375,      QE_R375,      lambda_(i), 'linear', 'extrap');
   QE_PMT_NEXT_(i)  = interp1(lambdaQE_PMT_NEXT,  QE_PMT_NEXT,  lambda_(i), 'linear', 'extrap');
   QE_CCD_(i)       = interp1(lambdaQE_CCD,       QE_CCD,       lambda_(i), 'linear', 'extrap');
   QE_H7826_(i)     = interp1(lambdaQE_H7826,     QE_H7826,     lambda_(i), 'linear', 'extrap');
   QE_H7826_01_(i)  = interp1(lambdaQE_H7826_01,  QE_H7826_01,  lambda_(i), 'linear', 'extrap');
   QE_H10425_(i)    = interp1(lambdaQE_H10425,    QE_H10425,    lambda_(i), 'linear', 'extrap');
   QE_H10425_01_(i) = interp1(lambdaQE_H10425_01, QE_H10425_01, lambda_(i), 'linear', 'extrap');
   QE_H10426_(i)    = interp1(lambdaQE_H10426,    QE_H10426,    lambda_(i), 'linear', 'extrap');
   QE_H10426_01_(i) = interp1(lambdaQE_H10426_01, QE_H10426_01, lambda_(i), 'linear', 'extrap');   
   QE_R1926_(i)     = interp1(lambdaQE_R1926,     QE_R1926,     lambda_(i), 'linear', 'extrap');
   QE_R7378_(i)     = interp1(lambdaQE_R7378,     QE_R7378,     lambda_(i), 'linear', 'extrap');   
   Transm_(i)       = interp1(lambdaQE_Transm,    Transm,       lambda_(i), 'linear', 'extrap');   
   %temp
   Yield100CF4_(i)  = interp1(lambda_CF4_tmp,     Yield100CF4,  lambda_(i), 'linear', 0);
end
% 
QE_H7826_     (QE_H7826_<0)     = 0;
QE_H7422_40_  (QE_H7422_40_<0)  = 0;
QE_H10721_01_ (QE_H10721_01_<0) = 0;
QE_H10721_20_ (QE_H10721_20_<0) = 0;
QE_H11461_01_ (QE_H11461_03_<0) = 0;
QE_R375_      (QE_R375_<0)      = 0;
QE_PMT_NEXT_  (QE_PMT_NEXT_<0)  = 0;
QE_CCD_       (QE_CCD_<0)       = 0;
QE_H7826_01_  (QE_H7826_01_<0)  = 0;
QE_H10425_    (QE_H10425_<0)    = 0;
QE_H10425_01_ (QE_H10425_01_<0) = 0;
QE_H10426_    (QE_H10426_<0)    = 0;
QE_H10426_01_ (QE_H10426_01_<0) = 0;
QE_R1926_     (QE_R1926_<0)     = 0;
QE_R7378_     (QE_R7378_<0)     = 0;
Transm_       (Transm_<0)       = 0;
Transm_       (lambda_<lambdaCut) = 0;

figure;
plot(lambda_, Transm_, '-');

Yield67Ar_  = Yield67Ar_/sum(Yield67Ar_);
Yield60He_  = Yield60He_/sum(Yield60He_);
Yield100CF4_ = Yield100CF4_/sum(Yield100CF4_);

Yield60Ar_  = Yield60Ar_/sum(Yield60Ar_);
Yield70Ar_  = Yield70Ar_/sum(Yield70Ar_);
Yield80Ar_  = Yield80Ar_/sum(Yield80Ar_);
Yield90Ar_  = Yield90Ar_/sum(Yield90Ar_);
Yield95Ar_  = Yield95Ar_/sum(Yield95Ar_);
Yield98Ar_  = Yield98Ar_/sum(Yield98Ar_);

figure; 
subplot(2,2,1); title('Light Yield'); hold on;
plot(lambda_, Yield67Ar_,   '-');
plot(lambda_, Yield95Ar_,   'r-');
plot(lambda_, Yield60He_, 'g-'); box;
legend('Ar/CF4(67/33)', 'Ar/CF4(95/5))', 'He/CF4(60/40)');

subplot(2,2,2); title('QE'); hold on;
plot(lambdaQE_H7826,     QE_H7826,     'r-');
plot(lambdaQE_H7826_01,  QE_H7826_01,  'r--');
plot(lambdaQE_R1926,     QE_R1926,     'b-');
plot(lambdaQE_R7378,     QE_R7378,     'b--');
plot(lambdaQE_H10426,    QE_H10426,    'g-');
plot(lambdaQE_H10426_01, QE_H10426_01, 'g--');
plot(lambdaQE_PMT_NEXT,  QE_PMT_NEXT,  'k-');
plot(lambdaQE_R375,      QE_R375,      'c')
plot(lambdaQE_CCD,       QE_CCD,       'k--'); box;
legend('H7826', 'H7826-01', 'R1926', 'R7378', 'H10426', 'H10426-01', 'NEXT', 'R375', 'CCD');

subplot(2,2,3); title('T*Y*QE Ar/CF4(67/33)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_,     'r-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_01_,  'r--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R1926_,     'b-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R7378_,     'b--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_,    'g-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_01_, 'g--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_,  'k-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R375_,      'c-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_CCD_,       'k--'); box;

subplot(2,2,4); title('T*Y*QE*Omega(10cm) Ar/CF4(67/33)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_     .* pi*(1.5/2)^2/(4*pi*10^2), 'r-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_01_  .* pi*(1.5/2)^2/(4*pi*10^2), 'r--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R1926_     .* pi*(2.2/2)^2/(4*pi*10^2), 'b-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R7378_     .* pi*(2.2/2)^2/(4*pi*10^2), 'b--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_    .* pi*(2.5/2)^2/(4*pi*10^2), 'g-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_01_ .* pi*(2.5/2)^2/(4*pi*10^2), 'g--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_  .* pi*(7.5/2)^2/(4*pi*10^2), 'k-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R375_      .* pi*(7.5/2)^2/(4*pi*10^2), 'c-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_CCD_       .* 1/16/(10+1)^2,            'k--'); box;

QE_67Ar_R1926   = sum(Yield67Ar_   .* Transm_ .* QE_R1926_)
QE_95Ar_R1926   = sum(Yield95Ar_   .* Transm_ .* QE_R1926_)
QE_60He_R1926   = sum(Yield60He_   .* Transm_ .* QE_R1926_)
QE_100CF4_R1926 = sum(Yield100CF4_ .* Transm_ .* QE_R1926_)

QE_67Ar_R7378   = sum(Yield67Ar_   .* Transm_ .* QE_R7378_)
QE_95Ar_R7378   = sum(Yield95Ar_   .* Transm_ .* QE_R7378_)
QE_60He_R7378   = sum(Yield60He_   .* Transm_ .* QE_R7378_)
QE_100CF4_R7378 = sum(Yield100CF4_ .* Transm_ .* QE_R7378_)

figure; title('yields'); hold on; 
plot(lambda_, Yield67Ar_  , '-');
plot(lambda_, Yield95Ar_ , 'r-');
plot(lambda_, Yield60He_ , 'g-');
plot(lambda_, Yield100CF4_ , 'k-'); box;
legend('Ar/CF4 67/33', 'Ar/CF4 95/5', 'He/CF4 0/60', 'CF4 100');

figure; title('T*Y*QE (R1926)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R1926_ , '-');
plot(lambda_, Yield95Ar_ .* Transm_ .* QE_R1926_ , 'r-');
plot(lambda_, Yield60He_ .* Transm_ .* QE_R1926_ , 'g-');
plot(lambda_, Yield100CF4_ .* Transm_ .* QE_R1926_ , 'k-'); box;
legend('Ar/CF4 67/33', 'Ar/CF4 95/5', 'He/CF4 0/60', 'CF4 100');

figure; title('T*Y*QE (R7378)'); hold on; 
plot(lambda_, Yield80Ar_ .* Transm_ .* QE_R7378_ , '-');
plot(lambda_, Yield95Ar_ .* Transm_ .* QE_R7378_ , 'r-');
plot(lambda_, Yield60He_ .* Transm_ .* QE_R7378_ , 'g-');
plot(lambda_, Yield100CF4_ .* Transm_ .* QE_R7378_ , 'k-'); box;
legend('Ar/CF4 80/20', 'Ar/CF4 95/5', 'He/CF4 0/60', 'CF4 100');


mosaic;