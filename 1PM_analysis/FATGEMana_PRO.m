%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    script for analizing FAT-GEM data with 1 PMT (DGD, 15/05/2021)   %%%%%%%
%%%%%%%                 adapted from (DGD, 20/04/2019)                      %%%%%%%
%%%%%%%       (for analysis of 2019 data use script in FALCON/ANALYSIS/OLD) %%%%%%%
%%%%%%%                 last version (DGD, 10/12/2022)                      %%%%%%%
%%%%%%%                                                                     %%%%%%%
%%%%%%%     -Data triggering for FAT-GEMs might be difficult without        %%%%%%%
%%%%%%%     limiting bandwidth at the scope (triggering in single-photons   %%%%%%%
%%%%%%%     creates trouble). It can be overcome with bandwidth reduction/  %%%%%%%
%%%%%%%     smart triggering in the scope (nCh =2), or triggering directly  %%%%%%%
%%%%%%%     on signal and later filtering events.                           %%%%%%%
%%%%%%%     -Smoothening and wvf-alignment implemented for easier analysis. %%%%%%%
%%%%%%%     -Event viewer implemented at the end of the script. Simply type %%%%%%%
%%%%%%%     the event indext to start, execute the entire block and then    %%%%%%%
%%%%%%%     keep pressing any button to scroll over events.                 %%%%%%%
%%%%%%%                                                                     %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;

% Path to aux scripts:

addpath [loc_path, 'HOME_RareEventsGroup\DiegoR\BasicScripts']
addpath [loc_path, 'HOME_RareEventsGroup\DiegoR\BasicScripts\_COMMON_SCRIPTS']

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% I. INITIALIZATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DataType           = 'CAEN';                                                %'CAEN' // 'OSCI' [data taken with CAEN card or with Tektronix scope]
Nsave              = 50000;                                                 %Number of events saved in workspace when wave-files are first read [required by ReadDataCAEN_PRO]

FILE               = 'wave';
nCh                = 1;                                                     %Indicate number of channels recorded                             [choice stored and might be overridden by SetRun]
isTriggerSignalIn1 = 0;                                                     %Indicate if trigger signal is, allows for smart external trigger [choice stored and might be overridden by SetRun]
FATGEMsetrun_PRO;                                                           %Set run and basic parameters  [might override previous choices]

minYForWVF  = -20;                                                          %Lowest value for wvf representation
Arange      = -10:0.5:500;                                                  %Range for amplitude representation (mV)
Erange      =  0:0.25:12;                                                   %Range for energy fit (spectrum)
Qrange      = -10:5:2000;                                                   %Range for charge representation (pC)
Wrange      = -0.1:0.01:2;                                                  %Range for width representation (us)
FATGEMini_PRO;                                                              %Set range for representation range and define cuts [might override previous choices]

if(isempty(Qth)), Qth = 0; end
if(isempty(Ath)), Ath = 0; end 

NsigmasThres = 3;                                                           %Threshold for waveform in number of sigmas (nsigmas over noise level)

isTriggerSetAutomatically = 1;                                              %Find trigger parameters automatically (default)
trigTime                  = 3;                                              %Put trigger parameters by hand (only if isTriggerSetAutomatically = 0)
trigWindow                = [1.8750, 4.4125];                               %Put trigger parameters by hand (only if isTriggerSetAutomatically = 0)

DoConvANDalign            = 1;                                              %reduce bandwidth and align waveforms for simpler analysis

isSmearedA    = 1*0;                                                        %smears amplitude
isSmearedT    = 1*0;                                                        %smears time

Rin           = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ib. READOUT (this part is automatic (no need to touch))             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(DataType,'CAEN')), ReadDataCAEN_PRO; end

if(isTriggerSignalIn1 == 1)
    data_       = data;        
    data        = [];
    data(1,:,:) = data_(2,:,:);    
    dataSize    = size(data);
end

Nchannels = dataSize(1);                                                    %number of PMTs
Nbins     = dataSize(2);                                                    %number of time channels
Nevts     = dataSize(3);                                                    %number of events

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% II. START ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%II.a) DEFINE WVF WINDOW AND TRIGGER TIME IF NOT AUTOMATIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AmpMeanWVF = zeros(1,Nchannels); TimeMeanWVF = zeros(1,Nchannels); tMinWVF = zeros(1,Nchannels); tMaxWVF = zeros(1,Nchannels); %#ok<PREALL>
 
data = squeeze(data(1,:,:));

if(isTriggerSetAutomatically)                                               %otherwise chose values set by hand in ini file
    MeanWVF             = mean(data');                                      %#ok<UDIM>
    PedMeanWVF          = MeanWVF(1);
    MeanWVF             = MeanWVF - PedMeanWVF;
    AmpMeanWVF          = max(MeanWVF);
    TimeMeanWVF         = time(MeanWVF==AmpMeanWVF);
    
    [~, ~, tMinWVF, tMaxWVF, ~, ~, ~, ~, ~, ~, ~, ~] = PSA_Nausicaa0(MeanWVF, time, AmpMeanWVF/10, TimeMeanWVF);
end
trigTime      = mean(TimeMeanWVF);
trigWindow(1) = 0.8*min(tMinWVF);
trigWindow(2) = 1.2*min(tMaxWVF);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%II.b) SUBTRACT PEDESTAL TO WVF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pedMEAN = zeros(Nchannels, Nevts); pedSTD = zeros(Nchannels, Nevts);      

data__  = zeros(size(data));                                                 %pedestal subtracted data

for i=1:Nevts
    dataTMP       = data(:,i);
    data_ped(:,i) = dataTMP(time<trigWindow(1));                  %#ok<SAGROW>
    
    pedMEAN(i)    = mean(data_ped(:,i));                          %used for pedestal subtraction
    pedSTD(i)     = std(data_ped(:,i));                           %used for threshold determination in PSA
    
    data__(:,i)   = data(:,i) - pedMEAN(i);                       %pedestal subtraction
end

if(DoConvANDalign)
    smoothSIGMA    = 200;
    smoothFunction = exp(-(time - 1000).^2/(2*smoothSIGMA));
    smoothFunction = 1/TBin*smoothFunction/(sum(smoothFunction));
    for i=1:Nevts
        if(round(i/100) == i/100), i=i; end
        dataTEMP = MTLconvol(data__(:,i), time, smoothFunction, time, 'tt');
        ti = mean(time(find(dataTEMP == max(dataTEMP))));
        tf = trigTime;
        dataTEMP = TimeAlign(time, dataTEMP, ti, tf);
        data__(:,i) = dataTEMP;
        
        data_ped      = dataTEMP(1:10);                                     %pedestal estimate   
        pedMEAN(i)    = mean(data_ped);                                     %used for pedestal subtraction
        pedSTD(i)     = std(data_ped);                                      %used for threshold determination in PSA
 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%II.d) PLOT WVFS FOR COARSE ASSESSMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure; hold on;
dataMean        = mean(squeeze(data__(:,:))');                   %#ok<UDIM>
plot(time, data__(:,1)); plot(time, data__(:,2), 'r'); plot(time, data__(:,3), 'g'); plot(time, data__(:,4), 'c'); plot(time, data__(:,5), 'k');
maxYForWVF = mean(max(data__(:,1:5)))*1.5; if(maxYForWVF<minYForWVF), maxYForWVF = -minYForWVF; end
plot(time, dataMean,'k'); lastline('LineWidth',3);
line([trigWindow(1), trigWindow(1)], [-500,500]); lastline('color','k');  lastline('LineStyle','--');
line([trigWindow(2), trigWindow(2)], [-500 500]); lastline('color','k');  lastline('LineStyle','--');
line([trigTime,      trigTime],      [-500 500]); lastline('color','r');
title('PMT wvfs');
yaxis(-2, maxYForWVF); xlabel('time[ns]'); ylabel('A[mV]'); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% III. 'DIGITIZATION/PULSE SHAPE ANALYSIS': %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%NOTE: threshold set at about x3 the average sigma (this should be tunable)

Qevt    = zeros(1, Nevts); Aevt    = zeros(1, Nevts); Wevt    = zeros(1, Nevts); Tevt = zeros(1, Nevts); 
Baseevt = zeros(1, Nevts); Trevt   = zeros(1, Nevts); Tfevt   = zeros(1, Nevts);
Assyevt = zeros(1, Nevts); stdTevt  = zeros(1, Nevts);

Vth = NsigmasThres*mean(pedSTD);

for i=1:Nevts
    F = data__(:,i)';  t = time;
    [Qevt(i), Aevt(i), ~, ~, Wevt(i), Tevt(i), ~, ~, Baseevt(i), Trevt(i), Tfevt(i), ~, Assyevt(i), stdTevt(i)] = PSA_Nausicaa0(F, t, Vth, trigTime);
    
    if(isSmearedA), Aevt(i)    = Aevt(i)    - ABin/2 + ABin*rand; end       %#ok<BDLGI>
    if(isSmearedT)                                                          %#ok<BDLGI>
        Wevt(i)  = Wevt(i)  - TBin/2 + TBin*rand;
        Trevt(i) = Trevt(i) - TBin/2 + TBin*rand;
        Tfevt(i) = Tfevt(i) - TBin/2 + TBin*rand;
    end
    
    Qevt(i)    = Qevt(i) * TBin / Rin;                   %[mV x ns = pC]
    
end

Wevt = Wevt/1000; Tevt = Tevt/1000; Trevt = Trevt/1000; Tfevt = Tfevt/1000; stdTevt = stdTevt/1000; %convert to us
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% IV. 'CUTS': %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Icut_SOFT = stdTevt>stdMIN   &   abs(Assyevt) < AssyevtMAX  &  Baseevt      < BasevtMAX      &  Aevt< AMAX;
Icut_HARD = stdTevt>stdMIN   &   stdTevt      < stdMAX   &  abs(Assyevt) < AssyevtMAX_HC  &  Baseevt < BasevtMAX_HC   &   Aevt< AMAX   &   Tevt < TevtMAX  &  Tevt > TevtMIN  & ...
    Trevt > TrevtMIN & Tfevt > TfevtMIN & Qevt > Qth & Aevt > Ath;

%PLOT event variables for different cuts
figure;

subplot(2,4,1); hold on;
hist1D(Aevt,             -10:1:500);
hist1D(Aevt(Icut_SOFT),  -10:1:500); lastline('color','r'); lastline('LineStyle','--');
hist1D(Aevt(Icut_HARD),  -10:1:500); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Aevt,             -10:1:500);
line([AMAX, AMAX],[0,max(N)]); lastline('color','r'); lastline('LineStyle','--');
line([AMAX, AMAX],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('amplitude [mV]'); ylabel('entries'); title ('amplitude'); box on;

subplot(2,4,2); hold on;
hist1D(Qevt,             Qrange);
hist1D(Qevt(Icut_SOFT),  Qrange); lastline('color','r'); lastline('LineStyle','--');
hist1D(Qevt(Icut_HARD),  Qrange); lastline('color','g'); lastline('LineStyle',':');
xlabel('charge [pC]'); ylabel('entries'); title ('charge'); box on;

subplot(2,4,3); hold on;
hist1D(stdTevt,             -0.02:0.01:2);
hist1D(stdTevt(Icut_SOFT),  -0.02:0.01:2); lastline('color','r'); lastline('LineStyle','--');
hist1D(stdTevt(Icut_HARD),  -0.02:0.01:2); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(stdTevt,           -0.02:0.01:2);
line([stdMIN, stdMIN],[0,max(N)]); lastline('color','r'); lastline('LineStyle','--');
line([stdMAX, stdMAX],[0,max(N)]); lastline('color','r'); lastline('LineStyle','--');
line([stdMIN, stdMIN],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
line([stdMAX, stdMAX],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('width [us]'); ylabel('entries'); title ('width (std. dev.)'); box on;

subplot(2,4,4); hold on;
hist1D(Tevt,               0:0.01:12);
hist1D(Tevt(Icut_SOFT),    0:0.01:12); lastline('color','r'); lastline('LineStyle','--');
hist1D(Tevt(Icut_HARD),    0:0.01:12); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Tevt,             0:0.01:12);
line([TevtMIN, TevtMIN],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
line([TevtMAX, TevtMAX],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
legend('no cuts', 'soft cuts', 'hard cuts');
xlabel('mean time [us]'); ylabel('entries'); title ('mean time'); box on;

subplot(2,4,5); hold on;
hist1D(Trevt,             -0.01:0.01:2);
hist1D(Trevt(Icut_SOFT),  -0.01:0.01:2); lastline('color','r'); lastline('LineStyle','--');
hist1D(Trevt(Icut_HARD),  -0.01:0.01:2); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Trevt,           -0.01:0.01:2);
line([TrevtMIN, TrevtMIN],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('rise time [us]'); ylabel('entries'); title ('rise time'); box on;

subplot(2,4,6); hold on;
hist1D(Tfevt,             -0.01:0.01:2);
hist1D(Tfevt(Icut_SOFT),  -0.01:0.01:2); lastline('color','r'); lastline('LineStyle','--');
hist1D(Tfevt(Icut_HARD),  -0.01:0.01:2); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Tfevt,           -0.01:0.01:2);
line([TfevtMIN, TfevtMIN],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('fall time [us]'); ylabel('entries'); title ('fall time'); box on;

subplot(2,4,7); hold on;
hist1D(Baseevt,             0:0.01:4);
hist1D(Baseevt(Icut_SOFT),  0:0.01:4); lastline('color','r'); lastline('LineStyle','--');
hist1D(Baseevt(Icut_HARD),  0:0.01:4); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Baseevt,           0:0.01:4);
line([BasevtMAX,      BasevtMAX],[0,max(N)]);    lastline('color','r'); lastline('LineStyle','--');
line([BasevtMAX_HC,   BasevtMAX_HC],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('base fluct [mV]'); ylabel('entries'); title ('base fluct'); box on;

subplot(2,4,8); hold on;
hist1D(Assyevt,             -10:0.1:10);
hist1D(Assyevt(Icut_SOFT),  -10:0.1:10); lastline('color','r'); lastline('LineStyle','--');
hist1D(Assyevt(Icut_HARD),  -10:0.1:10); lastline('color','g'); lastline('LineStyle',':');
N=hist1D(Assyevt,           -10:0.1:10);
line([AssyevtMAX,      AssyevtMAX],[0,max(N)]);    lastline('color','r'); lastline('LineStyle','--');
line([-AssyevtMAX,    -AssyevtMAX],[0,max(N)]);    lastline('color','r'); lastline('LineStyle','--');
line([AssyevtMAX_HC,   AssyevtMAX_HC],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
line([-AssyevtMAX_HC, -AssyevtMAX_HC],[0,max(N)]); lastline('color','g'); lastline('LineStyle',':');
xlabel('pulse assymetry [a.u.]'); ylabel('entries'); title ('pulse assym.'); box on;

%PLOT correlated variables
figure; hold on;
plot(stdTevt,           Qevt,'.');
plot(stdTevt(Icut_HARD),   Qevt(Icut_HARD),  'sg');
xaxis(min(Wrange), max(Wrange)); yaxis(min(Qrange), max(Qrange));
xlabel('width[us]'); ylabel('charge [pC]');
legend('no cuts', 'with cuts'); box on;
figure; hold on;
plot(Tfevt,          Qevt,'.');
plot(Tfevt(Icut_HARD),  Qevt(Icut_HARD),  'sg');
xaxis(min(Wrange), max(Wrange)); yaxis(min(Qrange), max(Qrange));
xlabel('tfall[us]'); ylabel('charge [pC]');
legend('no cuts', 'with cuts'); box on;
figure; hold on;
plot(Trevt,          Qevt,'.');
plot(Trevt(Icut_HARD),  Qevt(Icut_HARD),  'sg');
xaxis(min(Wrange), max(Wrange)); yaxis(min(Qrange), max(Qrange));
xlabel('trise[us]'); ylabel('charge [pC]');
legend('no cuts', 'with cuts'); box on;
figure; hold on;
plot(Tevt,           Qevt,'.');
plot(Tevt(Icut_HARD),   Qevt(Icut_HARD),  'sg');
yaxis(min(Qrange), max(Qrange));
xlabel('tevt[us]'); ylabel('charge [pC]');
legend('no cuts', 'with cuts'); box on;
figure; hold on;
plot(Aevt,           Qevt,'.');
plot(Aevt(Icut_HARD),   Qevt(Icut_HARD),  'sg');
xaxis(min(Arange), max(Arange)); yaxis(min(Qrange), max(Qrange));
xlabel('amplitude [mV]'); ylabel('charge [pC]');
legend('no cuts', 'with cuts'); box on;

%PLOT transient behaviour
figure; hold on;
plot(Qevt(Icut_HARD),'.');  yaxis(min(Qrange),max(Qrange)); title('drift'); xlabel('event number'); ylabel('event charge [pC]'); box;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% V. 'UNCORRELATE': %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Here tentatively uncorrelated with amplitude, although this is clearly not a valid approach
%
Ievt = 1:length(Qevt(Icut_HARD));
[Qevt_u, corr] = uncorr_D(Qevt(Icut_HARD), Ievt, 2, 't');
Qevt_u = Qevt_u + mean(Qevt(Icut_HARD));
plot(Qevt_u, 'or');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% VI. 'FIT AND SUMMARY': %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% SOFT CUTS %%%%%%%%

[N,X] = hist1D(Qevt(Icut_SOFT),  Qrange);
QFe55 = X(min(find(N==max(N))));                                            %#ok<MXFND>
FCal  = 5.9/QFe55;                                                          %Calibrate distributions based on cut one

%Fit distribution before cuts
QevtCal   = Qevt * FCal;

[Ndata, Qdata] = hist1D(QevtCal, Erange);                                
sigmaNdata     = sqrt(Ndata) + 0.0001;                                      %#ok<*NASGU>

par0(1)   = max(Ndata)/2;   %Amp
par0(2)   = 6;              %Mean 
par0(3)   = 5;              %Sigma 

ub(1)     = 2*max(Ndata);
ub(2)     = 6.5;
ub(3)     = 15;
lb(1)     = 0;
lb(2)     = 5.5;
lb(3)     = 0.2;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
par = lsqnonlin(@GFitterNausicaaFe55,par0,lb,ub,options);

figure; 
subplot(2,3,1); hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Fe_in_Xe_spectrum(par, Erange, 'Energy resolution'); title('raw');
xlabel('energy [keV]');
ResRaw     = par(3)*2.35/par(2);
QatPeakRaw = par(2)/FCal;

%Fit distribution after soft cuts
QevtCal   = Qevt(Icut_SOFT) * FCal;
[Ndata, Qdata] = hist1D(QevtCal, Erange);                                
sigmaNdata     = sqrt(Ndata) + 0.0001;                                      

par0(1)   = max(Ndata)/2;   %Amp
par0(2)   = 6;              %Mean 
par0(3)   = 2;              %Sigma 

ub(1)     = 2*max(Ndata);
ub(2)     = 10;
ub(3)     = 15;
lb(1)     = 0;
lb(2)     = 4;
lb(3)     = 0.2;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
par = lsqnonlin(@GFitterNausicaaFe55,par0,lb,ub,options);

subplot(2,3,2); hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Fe_in_Xe_spectrum(par, Erange, 'Energy resolution'); title('soft cut');
xlabel('energy [keV]');

ResCut     = par(3)*2.35/par(2);
QatPeakCut = par(2)/FCal;

%Fit distribution after soft cuts and uncorrelate
QevtCal   = Qevt_u * FCal;
[Ndata, Qdata] = hist1D(QevtCal, Erange);                                
sigmaNdata     = sqrt(Ndata) + 0.0001;                                      

par0(1)   = max(Ndata)/2;   %Amp
par0(2)   = 6;              %Mean 
par0(3)   = 2;              %Sigma 

ub(1)     = 2*max(Ndata);
ub(2)     = 10;
ub(3)     = 15;
lb(1)     = 0;
lb(2)     = 4;
lb(3)     = 0.2;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
par = lsqnonlin(@GFitterNausicaaFe55,par0,lb,ub,options);

subplot(2,3,3); hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Fe_in_Xe_spectrum(par, Erange, 'Energy resolution'); title('soft cut and uncorr');
xlabel('energy [keV]');

ResCutUnCorr     = par(3)*2.35/par(2);
QatPeakCutUnCorr = par(2)/FCal;

%%%%%%%% HARD CUTS %%%%%%%%

[N,X] = hist1D(Qevt(Icut_HARD),  Qrange);                                   %Calibrate distributions based on cut one
QFe55 = X(min(find(N==max(N))));                                            %#ok<MXFND>
FCal  = 5.9/QFe55;

%Fit distribution after hard cuts
QevtCal   = Qevt(Icut_HARD) * FCal;
[Ndata, Qdata] = hist1D(QevtCal, Erange);                                
sigmaNdata     = sqrt(Ndata) + 0.0001;                                      

par0(1)   = max(Ndata)/2;   %Amp
par0(2)   = 6;              %Mean 
par0(3)   = 2;              %Sigma 

ub(1)     = 2*max(Ndata);
ub(2)     = 10;
ub(3)     = 15;
lb(1)     = 0;
lb(2)     = 4;
lb(3)     = 0.2;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
par = lsqnonlin(@GFitterNausicaaFe55,par0,lb,ub,options);

subplot(2,3,5); hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Fe_in_Xe_spectrum(par, Erange, 'Energy resolution'); title('hard cut');
xlabel('energy [keV]');

ResHCut     = par(3)*2.35/par(2);
QatPeakHCut = par(2)/FCal;

%Fit distribution after hard cuts and uncorrelate
Ievt = 1:length(Qevt(Icut_HARD));
[Qevt_uHCut, corr] = uncorr_D(Qevt(Icut_HARD), Ievt, 2);
Qevt_uHCut = Qevt_uHCut + mean(Qevt(Icut_HARD));
QevtCal            = Qevt_uHCut * FCal;
[Ndata, Qdata]     = hist1D(QevtCal, Erange);                                
sigmaNdata         = sqrt(Ndata) + 0.0001;                                      

par0(1)   = max(Ndata)/2;   %Amp
par0(2)   = 6;              %Mean 
par0(3)   = 2;              %Sigma 

ub(1)     = 2*max(Ndata);
ub(2)     = 10;
ub(3)     = 15;
lb(1)     = 0;
lb(2)     = 4;
lb(3)     = 0.2;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
par = lsqnonlin(@GFitterNausicaaFe55,par0,lb,ub,options);

subplot(2,3,6); hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Fe_in_Xe_spectrum(par, Erange, 'Energy resolution'); title('hard cut and uncorr');
xlabel('energy [keV]');

ResHCutUnCorr     = par(3)*2.35/par(2);
QatPeakHCutUnCorr = par(2)/FCal;


toc;
mosaic;
disp(' summary ');
disp(' ');
disp(['energy resolution (raw)              = ', num2str(ResRaw),        '    events: ', num2str(length(Qevt))]);
disp(['energy resolution (soft cut)         = ', num2str(ResCut),        '    events: ', num2str(length(Qevt(Icut_SOFT)))]);
disp(['energy resolution (soft cut + corr)  = ', num2str(ResCutUnCorr)]);
disp(['energy resolution (hard cut)         = ', num2str(ResHCut),       '    events: ', num2str(length(Qevt(Icut_HARD)))]);
disp(['energy resolution (hard cut + corr)  = ', num2str(ResHCutUnCorr)]);
disp(' ');
disp(['charge at peak [pC]                  = ', num2str(QatPeakHCut)]);
disp(' ');
disp(['mean width (raw)                     = ', num2str(mean(stdTevt))]);
disp(['mean width (soft cut)                = ', num2str(mean(stdTevt(Icut_SOFT)))]);
disp(['mean width (hard cut)                = ', num2str(mean(stdTevt(Icut_HARD)))]);
disp(' ');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% VII. 'PLOT AVERAGE WVF': %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vth = NsigmasThres*mean(pedSTD);

dataMEAN_all = zeros(1,Nbins);
dataMEAN_HC  = zeros(1,Nbins);                                              %average values under hard cuts
dataMEAN_nHC = zeros(1,Nbins);                                              %average values for the complementary of the hard cuts set

for i=1:Nevts
    F = data__(:,i)';  t = time;
    dataMEAN_all = dataMEAN_all + F;
    if(Icut_HARD(i)==1),     dataMEAN_HC  = dataMEAN_HC + F;
    elseif(Icut_HARD(i)==0), dataMEAN_nHC = dataMEAN_nHC + F;
    end  
end
dataMEAN_all  = dataMEAN_all/Nevts;
dataMEAN_HC   = dataMEAN_HC/length(find(Icut_HARD>0));                     
dataMEAN_nHC  = dataMEAN_nHC/length(find(Icut_HARD==0));                  

figure; hold on;
plot(time, dataMEAN_all, '.b');
plot(time, dataMEAN_HC,  '.r');
plot(time, dataMEAN_nHC, '.g');
legend('all', 'hard cuts', 'complementary of hard cuts');
xlabel('time[ns]'); box;
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% VIII. 'Event viewer': %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vth = NsigmasThres*mean(pedSTD); 

for i=1:Nevts
    i=i
    F = data__(:,i)';  t = time;
    PSA_Nausicaa0(F, t, Vth, trigTime, 'hey');
    Icut_SOFT(i)
    Icut_HARD(i)
    Qevt(i)
    Aevt(i)
    pause;
    close;
end