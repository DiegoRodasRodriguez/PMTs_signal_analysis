close all;
clear all;

%NOTE:
%Correct for effective gain?. Indeed the Portuguese group collects e-
%at the bottom, too. It uses a negatively biased grid.
%Uncertainty in solid angle?.

lambdaCut   = 350;
npr         = 300;
dGem2Plexi  = 0.2;%*0.3/0.2;    %<0.3
tPlexi      = 0.3;%*0.5/0.3;
dPlexiMesh  = 0.25;%*0.5/0.25;
dMeshPMT    = 0.25;%*0.5/0.25;
dPlexiUp2PMT = 1.9;

GainGEM     = 10000;
GainLastGEM = 2*(GainGEM)^0.3; %approx
npe         = 80;

dir  = 'E:\HOME\+PROJECTS\++++2015_2016_ATCERN\1aOpticalTPC\AnalysisAndrea\QE\';

lambda_ = 200:1:900;

Yield67Ar_  = LightSpectrumInArCF4(0.33, lambda_);
Yield60Ar_  = LightSpectrumInArCF4(0.4, lambda_);
Yield70Ar_  = LightSpectrumInArCF4(0.3, lambda_);
Yield80Ar_  = LightSpectrumInArCF4(0.2, lambda_);
Yield90Ar_  = LightSpectrumInArCF4(0.1, lambda_);
Yield95Ar_  = LightSpectrumInArCF4(0.05, lambda_);
Yield98Ar_  = LightSpectrumInArCF4(0.02, lambda_);
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

file = 'TransmPlexiGlass';
eval(['load ', dir, file]);
lambdaQE_Transm = U(:,1); Transm  = U(:,2)/100 -0.05;

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
   Transm_(i)       = interp1(lambdaQE_Transm,    Transm,       lambda_(i), 'linear', 'extrap');
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
Transm_       (Transm_<0)       = 0;
Transm_       (lambda_<330)     = 0;

figure;
plot(lambda_, Transm_, '-');

Yield67Ar_  = Yield67Ar_/sum(Yield67Ar_(lambda_>400));
Yield60He_  = Yield60He_/sum(Yield60He_(lambda_>400));

Yield60Ar_  = Yield60Ar_/sum(Yield60Ar_(lambda_>400));
Yield70Ar_  = Yield70Ar_/sum(Yield70Ar_(lambda_>400));
Yield80Ar_  = Yield80Ar_/sum(Yield80Ar_(lambda_>400));
Yield90Ar_  = Yield90Ar_/sum(Yield90Ar_(lambda_>400));
Yield95Ar_  = Yield95Ar_/sum(Yield95Ar_(lambda_>400));
Yield98Ar_  = Yield98Ar_/sum(Yield98Ar_(lambda_>400));

figure; 
subplot(2,2,1); title('Light Yield'); hold on;
plot(lambda_, Yield67Ar_,   '-');
plot(lambda_, Yield95Ar_,   'r-');
plot(lambda_, Yield60He_, 'g-'); box;
legend('Ar/CF4(67/33)', 'Ar/CF4(95/5))', 'He/CF4(60/40)');

subplot(2,2,2); title('QE'); hold on;
plot(lambdaQE_H7826,     QE_H7826,     'r-');
plot(lambdaQE_H7826_01,  QE_H7826_01,  'r--');
plot(lambdaQE_H10425,    QE_H10425,    'b-');
plot(lambdaQE_H10425_01, QE_H10425_01, 'b--');
plot(lambdaQE_H10426,    QE_H10426,    'g-');
plot(lambdaQE_H10426_01, QE_H10426_01, 'g--');
plot(lambdaQE_PMT_NEXT,  QE_PMT_NEXT,  'k-');
plot(lambdaQE_R375,      QE_R375,      'c')
plot(lambdaQE_CCD,       QE_CCD,       'k--'); box;
legend('H7826', 'H7826-01', 'H10425', 'H10425-01', 'H10426', 'H10426-01', 'NEXT', 'R375', 'CCD');

subplot(2,2,3); title('T*Y*QE Ar/CF4(67/33)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_,     'r-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_01_,  'r--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10425_,    'b-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10425_01_, 'b--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_,    'g-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_01_, 'g--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_,  'k-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R375_,      'c-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_CCD_,       'k--'); box;

subplot(2,2,4); title('T*Y*QE*Omega(10cm) Ar/CF4(67/33)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_     .* pi*(1.5/2)^2/(4*pi*10^2), 'r-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_01_  .* pi*(1.5/2)^2/(4*pi*10^2), 'r--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10425_    .* pi*(2.2/2)^2/(4*pi*10^2), 'b-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10425_01_ .* pi*(2.2/2)^2/(4*pi*10^2), 'b--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_    .* pi*(2.5/2)^2/(4*pi*10^2), 'g-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H10426_01_ .* pi*(2.5/2)^2/(4*pi*10^2), 'g--');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_  .* pi*(7.5/2)^2/(4*pi*10^2), 'k-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_R375_      .* pi*(7.5/2)^2/(4*pi*10^2), 'c-');
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_CCD_       .* 1/16/(10+1)^2,            'k--'); box;

QE_67Ar_H7826 = sum(Yield67Ar_ .* Transm_ .* QE_H7826_)
QE_95Ar_H7826 = sum(Yield95Ar_ .* Transm_ .* QE_H7826_)
QE_60He_H7826 = sum(Yield60He_ .* Transm_ .* QE_H7826_)

% QE_67Ar_H7826      = sum(Yield67Ar_ .* Transm_ .* QE_H7826_)
% QE_67Ar_H7826_01   = sum(Yield67Ar_ .* Transm_ .* QE_H7826_01_)
% QE_67Ar_H10425     = sum(Yield67Ar_ .* Transm_ .* QE_H10425_)
% QE_67Ar_H10425_01  = sum(Yield67Ar_ .* Transm_ .* QE_H10425_01_)
% QE_67Ar_H10426     = sum(Yield67Ar_ .* Transm_ .* QE_H10426_)
% QE_67Ar_H10426_01  = sum(Yield67Ar_ .* Transm_ .* QE_H10426_01_)
% QE_67Ar_PMT_NEXT   = sum(Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_)
% QE_67Ar_CCD        = sum(Yield67Ar_ .* Transm_ .* QE_CCD_)

QEomega_67Ar_H7826      = sum(Yield67Ar_ .* Transm_ .* QE_H7826_     .* pi*(1.5/2)^2/(4*pi*10^2))
QEomega_67Ar_H7826_01   = sum(Yield67Ar_ .* Transm_ .* QE_H7826_01_  .* pi*(1.5/2)^2/(4*pi*10^2))
QEomega_67Ar_H10425     = sum(Yield67Ar_ .* Transm_ .* QE_H10425_    .* pi*(2.2/2)^2/(4*pi*10^2))
QEomega_67Ar_H10425_01  = sum(Yield67Ar_ .* Transm_ .* QE_H10425_01_ .* pi*(2.2/2)^2/(4*pi*10^2))
QEomega_67Ar_H10426     = sum(Yield67Ar_ .* Transm_ .* QE_H10426_    .* pi*(2.5/2)^2/(4*pi*10^2))
QEomega_67Ar_H10426_01  = sum(Yield67Ar_ .* Transm_ .* QE_H10426_01_ .* pi*(2.5/2)^2/(4*pi*10^2))
QEomega_67Ar_PMT_NEXT   = sum(Yield67Ar_ .* Transm_ .* QE_PMT_NEXT_  .* pi*(7.5/2)^2/(4*pi*10^2))
QEomega_67Ar_CCD        = sum(Yield67Ar_ .* Transm_ .* QE_CCD_       .* 1/16/(10+1)^2)


figure; title('T*Y*QE (H7286)'); hold on; 
plot(lambda_, Yield67Ar_ .* Transm_ .* QE_H7826_ , '-');
plot(lambda_, Yield95Ar_ .* Transm_ .* QE_H7826_ , 'r-');
plot(lambda_, Yield60He_ .* Transm_ .* QE_H7826_ , 'g-'); box;
legend('Ar/CF4 67/33', 'Ar/CF4 95/5', 'Ar/CF4 100/0');

figure; title('T*Y*QE (R325)'); hold on; 
plot(lambda_, Yield80Ar_ .* Transm_ .* QE_R375_ , '-');
plot(lambda_, Yield95Ar_ .* Transm_ .* QE_R375_ , 'r-');
plot(lambda_, Yield60He_ .* Transm_ .* QE_R375_ , 'g-'); box;
legend('Ar/CF4 80/20', 'Ar/CF4 95/5', 'Ar/CF4 100/0');


QE_67Ar_H7826 = sum(Yield67Ar_ .* Transm_ .* QE_H7826_)
QE_60He_H7826 = sum(Yield60He_ .* Transm_ .* QE_H7826_)

QE_60Ar_H7826 = sum(Yield60Ar_ .* Transm_ .* QE_H7826_)
QE_70Ar_H7826 = sum(Yield70Ar_ .* Transm_ .* QE_H7826_)
QE_80Ar_H7826 = sum(Yield80Ar_ .* Transm_ .* QE_H7826_)
QE_90Ar_H7826 = sum(Yield90Ar_ .* Transm_ .* QE_H7826_)
QE_95Ar_H7826 = sum(Yield95Ar_ .* Transm_ .* QE_H7826_)
QE_98Ar_H7826 = sum(Yield98Ar_ .* Transm_ .* QE_H7826_)

%Omega = pi*(1.5/2)^2/(4*pi*(dGem2Plexi+tPlexi+dPlexiMesh+dMeshPMT)^2)
Omega = pi*(1.5/2)^2/(4*pi*(dGem2Plexi+dPlexiUp2PMT)^2)

nphPerPrimary_67 = 1/(QE_67Ar_H7826 * Omega) / npr * npe / GainGEM
nphPerPrimary_95 = 1/(QE_95Ar_H7826 * Omega) / npr * npe / GainGEM
nphPerPrimary_60 = 1/(QE_60He_H7826 * Omega) / npr * npe / GainGEM

%               06/08/2015
%               67%  -  95%;
% nphPerPrim = 0.0045 - 0.0076 (200nm);
% nphPerPrim = 0.0086 - 0.0101 (400nm);
% nphPerPrim = 0.0179 - 0.0286 (300nm & assuming at PMT distance);
% x times 2 if considering effective gain (i.e., electrons lost to GEM).
%
% estimate nphPerPrim = 0.015 +- 0.01;

QE_60Ar_R375 = sum(Yield60Ar_ .* Transm_ .* QE_R375_)
QE_70Ar_R375 = sum(Yield70Ar_ .* Transm_ .* QE_R375_)
QE_80Ar_R375 = sum(Yield80Ar_ .* Transm_ .* QE_R375_)
QE_90Ar_R375 = sum(Yield90Ar_ .* Transm_ .* QE_R375_)
QE_95Ar_R375 = sum(Yield95Ar_ .* Transm_ .* QE_R375_)
QE_98Ar_R375 = sum(Yield98Ar_ .* Transm_ .* QE_R375_)

mosaic;