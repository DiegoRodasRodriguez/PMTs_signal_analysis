%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script for analyzing single-photon PMT-data from LED            %%%%%
%%%%%                                                                 %%%%%
%%%%%                   (DGD 31/Aug/2020)                             %%%%%
%%%%%               (final version 10/Dec/2022)                       %%%%%
%%%%%                                                                 %%%%%
%%%%%   USE NOTES:                                                    %%%%%
%%%%%                                                                 %%%%%
%%%%%   - File read and saved into workspace for easier reanalysis.   %%%%%
%%%%%   - The program works automatically except for the analysis     %%%%%
%%%%%   time range, that needs to be specified after a first run      %%%%%
%%%%%   (this could be done automatically... but this implementation  %%%%%
%%%%%   also forces the user to visually check that everything is ok).%%%%%
%%%%%   - Windows chosen for analysis and pedestal estimate are       %%%%%  
%%%%%   shown in wvf plots.                                           %%%%%    
%%%%%   - Check parameters in section Ia before running.              %%%%%
%%%%%   - It is recommended to store the configuration in the _ini    %%%%%
%%%%%   file to ease future analysis and keep a log.                  %%%%%
%%%%%                                                                 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;

addpath [loc_path, 'HOME_RareEventsGroup\DiegoR\BasicScripts']
addpath [loc_path, 'HOME_RareEventsGroup\DiegoR\BasicScripts\_COMMON_SCRIPTS']

global Ndata Qdata sigmaNdata;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ia. INITIALIZATION (adjust your desired values here)                %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NevtsMax        = 10000;                                                    %Number of events to be analyzed
Nsave           = 30000;                                                    %Number of events saved in workspace when wave files are first read (required by ReadDataCAEN_PRO)                                      

DataType        = 'CAEN';                                                   %'CAEN' // 'OSCI' [data taken with CAEN card or with Tektronix scope]
nCh             = 2;                                                        %Indicate number of channels recorded [choice stored and might be overridden by AnaSinglePhoton_ini]
isLEDsignal1    = 1;                                                        %Indicate if LED signal is in data    [choice stored and might be overridden by AnaSinglePhoton_ini]
iPM             = 2;                                                        %Indicate PM channel in case several are stored
iTrig           = 1;                                                        %Indicate trigger channel (by default containing the digital signal biasing the LED)
FILE            = 'wave';                                                   %File names                           [choice stored and might be overridden by AnaSinglePhoton_ini]
AmplitudeWindow = [700,780];                                                %Window to look for signal amplitude  [choice stored and might be overridden by AnaSinglePhoton_ini]
ChargeWindow    = [700,780];                                                %Window to look for signal charge     [choice stored and might be overridden by AnaSinglePhoton_ini]
SinglePhotonIni;                                                            %Set analysis parameters in ini file  [might override previous choices]

iSignal1 = 1;                                                               %signal index for visualization plot (by default, chose the first 5)
iSignal2 = 2;
iSignal3 = 3;
iSignal4 = 4;
iSignal5 = 5;

isRangeAutomatic = 1*0;                                                     %Automatic range for representation (otherwise use values below)                                          
Amin = -1; Amax = 50;
Qmin = -2; Qmax = 6;
QBin         = 0.05;                                                        %[pC] (default is 50 fC), ABin is kept to the one marked by hardware, after the reader function (e.g. ReadDataCAEN_PRO)
isSmeared    = 1;                                                           %Smears amplitude within each bin (for nicer representation)

isAFit = 1;                                                                 %Try fitting amplitude distribution, if data quality allows (charge data is always fit)
                                                        
isSaveToFile = 1*0;                                                         %Save signal to file so to port to other analysis software
DIR_OUTPUT   = DIR;                                                         %Output file to save signals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ib. READ OUT (from this point onwards no need to touch)             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(DataType,'OSCI')), ReadDataOsciSingle; end                        %Read from Tektronix scope
if(strcmp(DataType,'CAEN')), ReadDataCAEN_PRO; end                          %Read from CAEN card (requires defining Nsave)

sizedata          = size(data);
isTriggerSignalIn = (length(sizedata)>2);
if(NevtsMax && ~isTriggerSignalIn)
    Nevts = NevtsMax;    
    data_       = data; data = [];
    data(:,:)   = data_(:, 1:Nevts);
elseif(NevtsMax && isTriggerSignalIn)
    Nevts = NevtsMax; 
    data_       = data; data = [];
    data(:,:,:) = data_(:, :, 1:Nevts);
end

if(isTriggerSignalIn && isLEDsignal1 ==1)
    data_       = data; data = [];
    data(1,:,:) = data_(2,:,:);    
    data(2,:,:) = data_(1,:,:);
    iPM          = 1;
    iTrig        = 2; 
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IIa. START ANALYSIS (OBTAIN WVF IN WINDOWS)                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_AWindow = (time>AmplitudeWindow(1) & time<AmplitudeWindow(2));
I_QWindow = (time>ChargeWindow(1)    & time<ChargeWindow(2));

time_Awindow = time(I_AWindow);
time_Qwindow = time(I_QWindow);  

SPsignal_Awindow = zeros(length(time_Awindow), Nevts); 
SPsignal_AwindowPed = zeros(length(time_Awindow), Nevts); SPsignal_AwindowPed2 = zeros(length(time_Awindow), Nevts);
SPsignal_Qwindow = zeros(length(time_Qwindow), Nevts);
SPsignal_QwindowPed = zeros(length(time_Qwindow), Nevts); SPsignal_QwindowPed2 = zeros(length(time_Qwindow), Nevts);
Triggersignal_AwindowPed = zeros(length(time_Awindow), Nevts);

indexes_Awindow = find(I_AWindow);
indexes_Qwindow = find(I_QWindow);

if(~isTriggerSignalIn)
    SPsignal = data(:,:);
elseif(isTriggerSignalIn)
    SPsignal      = squeeze(data(iPM,:,:));
    Triggersignal = squeeze(data(iTrig,:,:));
end

iminPedA = 1; imaxPedA = (indexes_Awindow(length(indexes_Awindow)) - indexes_Awindow(1)) + 1;
iminPedA2 = imaxPedA + 1; imaxPedA2 = iminPedA2 + (imaxPedA -iminPedA);
    
iminPedQ = 1; imaxPedQ = (indexes_Qwindow(length(indexes_Qwindow)) - indexes_Qwindow(1)) + 1;
iminPedQ2 = imaxPedQ+1; imaxPedQ2 = iminPedQ2 + (imaxPedQ-iminPedQ);
    
for i=1:Nevts
    SPsignal_ = SPsignal(:,i);
    
    SPsignal_Awindow(:,i)     = SPsignal_(I_AWindow);
    SPsignal_AwindowPed(:,i)  = SPsignal_(iminPedA :imaxPedA );
    SPsignal_AwindowPed2(:,i) = SPsignal_(iminPedA2:imaxPedA2);                     %As pedestal is defined event-by-event, another reference is needed for estimating the pedestal spread
    
    SPsignal_Qwindow(:,i)     = SPsignal_(I_QWindow);
    SPsignal_QwindowPed(:,i)  = SPsignal_(iminPedQ :imaxPedQ );
    SPsignal_QwindowPed2(:,i) = SPsignal_(iminPedQ2:imaxPedQ2);                     %As pedestal is defined event-by-event, another reference is needed for estimating the pedestal spread
    
    if(isTriggerSignalIn)
        Triggersignal_AwindowPed(:,i) = Triggersignal(1:(indexes_Awindow(length(indexes_Awindow))-indexes_Awindow(1))+1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IIb. START ANALYSIS (SUBTRACT PEDESTAL TO WVF)                      %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Mean_APed    = zeros(1, Nevts); 
Mean_TrigPed = zeros(1, Nevts); 

SPsignal_Awindow_ = zeros(size(SPsignal_Awindow)); SPsignal_AwindowPed_ = zeros(size(SPsignal_AwindowPed));
SPsignal_Qwindow_ = zeros(size(SPsignal_Qwindow)); SPsignal_QwindowPed_ = zeros(size(SPsignal_Qwindow));
SPsignal_         = zeros(size(SPsignal));         Triggersignal_       = zeros(size(SPsignal)); 

for i=1:Nevts
    Mean_APed(:,i)            = mean(SPsignal_AwindowPed(:,i));
    SPsignal_(:,i)            = SPsignal(:,i)             - Mean_APed(i);  
    SPsignal_Awindow_(:,i)    = SPsignal_Awindow(:,i)     - Mean_APed(i);  
    SPsignal_AwindowPed_(:,i) = SPsignal_AwindowPed2(:,i) - Mean_APed(i); 
    SPsignal_Qwindow_(:,i)    = SPsignal_Qwindow(:,i)     - Mean_APed(i); 
    SPsignal_QwindowPed_(:,i) = SPsignal_QwindowPed2(:,i) - Mean_APed(i);  
          
    Mean_TrigPed(:,i)         = mean(Triggersignal_AwindowPed(:,i)); 
    
    if(isTriggerSignalIn)
        Triggersignal_(:,i)       = Triggersignal(:,i)        - Mean_TrigPed(i);
    end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IIc. START ANALYSIS (PLOT WVFS)                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SPsignalMean          = mean(SPsignal(:,:)');                                     %#ok<*UDIM>
SPsignalMean_         = mean(SPsignal_(:,:)');                                    %#ok<*UDIM>
if(isTriggerSignalIn), TriggersignalMean     = mean(Triggersignal(:,:)');  end    %#ok<*UDIM>
if(isTriggerSignalIn), TriggersignalMean_    = mean(Triggersignal_(:,:)'); end    %#ok<*UDIM>

minSignal   = min(min(SPsignal(:,1:20)));       maxSignal   = max(max(SPsignal(:,1:20)));
if(isTriggerSignalIn), minTSignal  = min(min(Triggersignal(:,1:20)));  maxTSignal  = max(max(Triggersignal(:,1:20))); end
minSignal_  = min(min(SPsignal_(:,1:20)));      maxSignal_  = max(max(SPsignal_(:,1:20)));
if(isTriggerSignalIn), minTSignal_ = min(min(Triggersignal_(:,1:20))); maxTSignal_ = max(max(Triggersignal_(:,1:20))); end

figure;                                                                     %pedestal-uncorrected
if(isTriggerSignalIn), subplot(2,1,1); end
hold on;
plot(time, SPsignal(:,iSignal1)); plot(time, SPsignal(:,iSignal2), 'r'); plot(time, SPsignal(:,iSignal3), 'g');plot(time, SPsignal(:,iSignal4), 'c'); plot(time, SPsignal(:,iSignal5), 'k');plot(time, SPsignalMean,'k'); lastline('LineWidth',3);
line([AmplitudeWindow(1), AmplitudeWindow(1)], [minSignal, maxSignal]); lastline('color','b');
line([AmplitudeWindow(2), AmplitudeWindow(2)], [minSignal, maxSignal]); lastline('color','b');
line([ChargeWindow(1),    ChargeWindow(1)],    [minSignal, maxSignal]); lastline('color','r'); lastline('LineStyle','--');
line([ChargeWindow(2),    ChargeWindow(2)],    [minSignal, maxSignal]); lastline('color','r'); lastline('LineStyle','--');
line([time(iminPedA),     time(iminPedA)],     [minSignal, maxSignal]); lastline('color','g');
line([time(iminPedA2),    time(iminPedA2)],    [minSignal, maxSignal]); lastline('color','g');
line([time(imaxPedA),     time(imaxPedA)],     [minSignal, maxSignal]); lastline('color','g');
line([time(imaxPedA2),    time(imaxPedA2)],    [minSignal, maxSignal]); lastline('color','g');
line([time(iminPedQ),     time(iminPedQ)],     [minSignal, maxSignal]); lastline('color','k'); lastline('LineStyle','--');
line([time(iminPedQ2),    time(iminPedQ2)],    [minSignal, maxSignal]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ),     time(imaxPedQ)],     [minSignal, maxSignal]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ2),    time(imaxPedQ2)],    [minSignal, maxSignal]); lastline('color','k'); lastline('LineStyle','--');
plot(time, SPsignalMean,'k'); lastline('LineWidth',3);  %replot for clarity
title('SP wvf (inverted)'); xlabel('time[ns]'); ylabel('A[mV]'); xaxis(time(1), time(length(time))); 
legend('signal 1', 'signal 2', 'signal 3', 'signal 4', 'signal 5', 'avg'); box;
if(isTriggerSignalIn) 
    subplot(2,1,2); hold on;
    plot(time, Triggersignal(:,iSignal1)); plot(time, Triggersignal(:,iSignal2), 'r'); plot(time, Triggersignal(:,iSignal3), 'g');plot(time, Triggersignal(:,iSignal4), 'c'); plot(time, Triggersignal(:,iSignal5), 'k');
    line([AmplitudeWindow(1), AmplitudeWindow(1)], [minTSignal, maxTSignal]); lastline('color','b');
    line([ChargeWindow(1),    ChargeWindow(1)],    [minTSignal, maxTSignal]); lastline('color','r'); lastline('LineStyle','--');
    line([AmplitudeWindow(2), AmplitudeWindow(2)], [minTSignal, maxTSignal]); lastline('color','b');
    line([ChargeWindow(2),    ChargeWindow(2)],    [minTSignal, maxTSignal]); lastline('color','r'); lastline('LineStyle','--');
    plot(time, TriggersignalMean,'k'); lastline('LineWidth',3);
    title('trigger wvf');
    xlabel('time[ns]'); ylabel('A[mV]'); xaxis(time(1), time(length(time))); box;
end

figure;                                                                     %pedestal-corrected
if(isTriggerSignalIn), subplot(2,1,1); end
hold on;
plot(time, SPsignal_(:,iSignal1)); plot(time, SPsignal_(:,iSignal2), 'r'); plot(time, SPsignal_(:,iSignal3), 'g');plot(time, SPsignal_(:,iSignal4), 'c'); plot(time, SPsignal_(:,iSignal5), 'k');plot(time, SPsignalMean_,'k'); lastline('LineWidth',3);
line([AmplitudeWindow(1), AmplitudeWindow(1)], [minSignal_, maxSignal_]); lastline('color','b');
line([AmplitudeWindow(2), AmplitudeWindow(2)], [minSignal_, maxSignal_]); lastline('color','b');
line([ChargeWindow(1),    ChargeWindow(1)],    [minSignal_, maxSignal_]); lastline('color','r'); lastline('LineStyle','--');
line([ChargeWindow(2),    ChargeWindow(2)],    [minSignal_, maxSignal_]); lastline('color','r'); lastline('LineStyle','--');
line([time(iminPedA),     time(iminPedA)],     [minSignal_, maxSignal_]); lastline('color','g');
line([time(iminPedA2),    time(iminPedA2)],    [minSignal_, maxSignal_]); lastline('color','g');
line([time(imaxPedA),     time(imaxPedA)],     [minSignal_, maxSignal_]); lastline('color','g');
line([time(imaxPedA2),    time(imaxPedA2)],    [minSignal_, maxSignal_]); lastline('color','g');
line([time(iminPedQ),     time(iminPedQ)],     [minSignal_, maxSignal_]); lastline('color','k'); lastline('LineStyle','--');
line([time(iminPedQ2),    time(iminPedQ2)],    [minSignal_, maxSignal_]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ),     time(imaxPedQ)],     [minSignal_, maxSignal_]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ2),    time(imaxPedQ2)],    [minSignal_, maxSignal_]); lastline('color','k'); lastline('LineStyle','--');
plot(time, SPsignalMean_,'k'); lastline('LineWidth',3);  %replot for clarity
title('SP wvf (inverted and pedestal corrected)'); xlabel('time[ns]'); ylabel('A[mV]'); xaxis(time(1), time(length(time)));
legend('signal 1', 'signal 2', 'signal 3', 'signal 4', 'signal 5', 'avg'); box;
if(isTriggerSignalIn) 
    subplot(2,1,2); hold on;
    plot(time, Triggersignal_(:,iSignal1)); plot(time, Triggersignal_(:,iSignal2), 'r'); plot(time, Triggersignal_(:,iSignal3), 'g');plot(time, Triggersignal_(:,iSignal4), 'c'); plot(time, Triggersignal_(:,iSignal5), 'k');
    line([AmplitudeWindow(1), AmplitudeWindow(1)], [minTSignal_, maxTSignal_]); lastline('color','b');
    line([ChargeWindow(1),    ChargeWindow(1)],    [minTSignal_, maxTSignal_]); lastline('color','r'); lastline('LineStyle','--');
    line([AmplitudeWindow(2), AmplitudeWindow(2)], [minTSignal_, maxTSignal_]); lastline('color','b');
    line([ChargeWindow(2),    ChargeWindow(2)],    [minTSignal_, maxTSignal_]); lastline('color','r'); lastline('LineStyle','--');
    plot(time, TriggersignalMean_,'k'); lastline('LineWidth',3);
    title('trigger wvf (pedestal corrected)');
    xlabel('time[ns]'); ylabel('A[mV]'); xaxis(time(1), time(length(time))); box;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% III. 'DIGITIZATION/PULSE SHAPE ANALYSIS':Obtain amplitude and charge %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AMPLITUDE
Amp    = max(SPsignal_Awindow_(:,:));
AmpPed = max(SPsignal_AwindowPed_(:,:));
if(isSmeared)
    Amp    = Amp    - ABin/2 + ABin*rand(size(Amp));
    AmpPed = AmpPed - ABin/2 + ABin*rand(size(AmpPed)); 
end
%CHARGE
Rin = 50;
Charge    = squeeze(sum(SPsignal_Qwindow_(:,:)))    * TBin / Rin; %pC
ChargePed = squeeze(sum(SPsignal_QwindowPed_(:,:))) * TBin / Rin; %pC

if(isRangeAutomatic)                                                        %#ok<BDLGI>
    Arange     = 0:ABin:max(Amp);
    Qrange     = min(Charge):QBin:max(Charge);
    ArangePed  = 0:ABin:max(AmpPed);
    QrangePed  = min(ChargePed):QBin:max(ChargePed);
else
    Arange    = Amin:ABin:Amax;
    Qrange    = Qmin:QBin:Qmax;
    ArangePed = Amin:ABin:Amax;
    QrangePed = Qmin:QBin:Qmax;
end

figure;
subplot(2,1,1);
hold on;
hist1D(Amp,    Arange);
hist1D(AmpPed, Arange); lastline('color', 'r')
legend('amplitude', 'noise(in the same window)'); xaxis(min(Arange),max(Arange));
N = hist1D(Amp, Arange); yaxis(0,1.2*max(N));
xlabel('A[mV]'); ylabel(['entries (', num2str(ABin), 'mV bin)']); title ('amplitude'); box;
subplot(2,1,2);
hold on;
hist1D(Charge,    Qrange);
hist1D(ChargePed, Qrange);lastline('color', 'r')
legend('Charge', 'noise(in the same window)'); xaxis(min(Qrange),max(Qrange));
N = hist1D(Charge, Qrange); yaxis(0,1.2*max(N));
xlabel('Q[pC]');  ylabel(['entries (', num2str(QBin), 'pC bin)']); title ('charge'); box;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IV. 'FIT'                                                           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Charge distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ndata, Qdata] = hist1D(Charge,  Qrange);

sigmaNdata     = sqrt(Ndata) + 0.0001;                                      %#ok<NASGU>

par0(1)   = 1;    %Average number of photons
par0(2)   = 1000; %Amplitude of single-photon peak
par0(3)   = 1;    %Width of single-photon peak
par0(4)   = 1;    %Mean of single-photon peak
par0(5)   = 0.2;  %Width of pedestal
par0(6)   = 0;    %Mean of pedestal

lb(1)     = 0;     ub(1)     = 5;
lb(2)     = 50;    ub(2)     = 3000;
lb(3)     = 0.1;   ub(3)     = 2;
lb(4)     = 0.2;   ub(4)     = 2;
lb(5)     = 0.05;  ub(5)     = 1.5;
lb(6)     = -0.5;  ub(6)     = 0.5;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parQ = lsqnonlin(@GFitterSP,par0,lb,ub,options);                            %lsqnonlin requieres Optimization Toolbox add-on installed

figure; hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
SP_spectrum(parQ, Qrange, 'Single-Photon distribution (Q)'); xaxis(min(Qrange),max(Qrange));
xlabel('Q[pC]');  ylabel(['entries (', num2str(QBin), ' pC bin)']); title ('Charge in window');

%Charge pedestal distribution
[Ndata, Qdata] = hist1D(ChargePed,  QrangePed);

sigmaNdata     = sqrt(Ndata) + 0.0001;                                     

par0(1)   = 0;    %Mean position of pedestal
par0(2)   = 1;    %Width of pedestal
par0(3)   = 1000; %Amplitude of pedestal

lb(1)     = -1;    ub(1)     = 1;
lb(2)     =  0.01; ub(2)     = 1.5;
lb(3)     = 10;    ub(3)     = 10000;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parQped = lsqnonlin(@GFitterSPped,par0,lb,ub,options);

figure; hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
SP_spectrumPed(parQped, QrangePed, 'Pedestal distribution (Q)'); xaxis(min(Qrange),max(Qrange));
xlabel('Q[pC]');  ylabel(['entries (', num2str(QBin), ' pC bin)']); title ('Charge in pedestal');

%%%%%%%%%%%%%%%%%%%%%% Amplitude distribution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Ndata, Qdata] = hist1D(Amp, Arange);

sigmaNdata     = sqrt(Ndata) + 0.0001;                                      %#ok<NASGU>

par0(1)   = 1;    %Average number of photons
par0(2)   = 300;  %Amplitude of single-photon peak
par0(3)   = 10;   %Width of single-photon peak
par0(4)   = 10;   %Mean of single-photon peak
par0(5)   = 1;    %Width of pedestal
par0(6)   = 3;    %Mean of pedestal

lb(1)     = 0;     ub(1)     = 5;
lb(2)     = 5;     ub(2)     = 3000;
lb(3)     = 3;     ub(3)     = 30;
lb(4)     = 3;     ub(4)     = 30;
lb(5)     = 0;     ub(5)     = 5;
lb(6)     = 0;     ub(6)     = 5;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parA = lsqnonlin(@GFitterSP,par0,lb,ub,options);

figure; hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
SP_spectrum(parA, Arange, 'Single-Photon distribution (A)'); xaxis(min(Arange),max(Arange));
xlabel('A[mV]');  ylabel(['entries (', num2str(ABin), ' mV bin)']); title ('Amplitude in window');

%Charge pedestal distribution
[Ndata, Qdata] = hist1D(AmpPed,  ArangePed);

sigmaNdata     = sqrt(Ndata) + 0.0001;                                     

par0(1)   = 3;    %Mean position of pedestal
par0(2)   = 1;    %Width of pedestal
par0(3)   = 1000; %Amplitude of pedestal

lb(1)     = 0;    ub(1)     = 5;
lb(2)     = 0;    ub(2)     = 5;
lb(3)     = 10;   ub(3)     = 10000;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parAped = lsqnonlin(@GFitterSPped,par0,lb,ub,options);

figure; hold on;
plot(Qdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
SP_spectrumPed(parAped, ArangePed, 'Pedestal distribution (A)'); xaxis(min(Arange),max(Arange));
xlabel('A[mV]');  ylabel(['entries (', num2str(ABin), ' mV bin)']); title ('Amplitude in pedestal');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% V. 'PLOT AFTER SUBTRACTING CHARGE-PEDESTALS'                        %%%
%%%  (Cut for detailed study of non-zero events)                        %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

j=0;
for i=1:length(Amp)
    if(Charge(i) > (parQ(6) + 2*parQ(5)))                                   %2sigmas above pedestal                              
        j=j+1;
        Amp_cut(j)         = Amp(i);                                        %#ok<SAGROW>
        Charge_cut(j)      = Charge(i);                                     %#ok<SAGROW>
        SPsignal_cut_(:,j) = SPsignal_(:,i);                                %#ok<SAGROW>
    end
end

figure;
subplot(2,1,1);
hold on;
hist1D(Amp,     Arange);
hist1D(Amp_cut, Arange); lastline('color', 'g');
legend('before cuts','after cuts'); xaxis(min(Arange),max(Arange));
N = hist1D(Amp, Arange); yaxis(0,1.2*max(N));
xlabel('A[mV]'); ylabel(['entries (', num2str(ABin), 'mV bin)']); title ('amplitude'); box;
subplot(2,1,2);
hold on;
hist1D(Charge,     Qrange);
hist1D(Charge_cut, Qrange); lastline('color', 'g');
legend('before cuts','after cuts'); xaxis(min(Qrange),max(Qrange));
N = hist1D(Charge, Qrange); yaxis(0,1.2*max(N));
xlabel('Q[pC]');  ylabel(['entries (', num2str(QBin), 'pC bin)']); title ('Charge'); box;

SPsignalMean_cut          = mean(SPsignal_cut_(:,:)');                      %#ok<*UDIM>
minSignalCut = min(SPsignalMean_cut)*1.1;  maxSignalCut = max(SPsignalMean_cut)*1.1;

figure;
hold on;
plot(time, SPsignalMean_,'r');
plot(time, SPsignalMean_cut,'k'); lastline('LineWidth',3);
line([AmplitudeWindow(1), AmplitudeWindow(1)],  [minSignalCut, minSignalCut]); lastline('color','b');
line([ChargeWindow(1),    ChargeWindow(1)],     [minSignalCut, maxSignalCut]); lastline('color','r'); lastline('LineStyle','--');
line([AmplitudeWindow(2), AmplitudeWindow(2)],  [minSignalCut, maxSignalCut]); lastline('color','b');
line([ChargeWindow(2),    ChargeWindow(2)],     [minSignalCut, maxSignalCut]); lastline('color','r'); lastline('LineStyle','--');
line([time(iminPedA),     time(iminPedA)],      [minSignalCut, maxSignalCut]); lastline('color','g');
line([time(iminPedA2),    time(iminPedA2)],     [minSignalCut, maxSignalCut]); lastline('color','g');
line([time(imaxPedA),     time(imaxPedA)],      [minSignalCut, maxSignalCut]); lastline('color','g');
line([time(imaxPedA2),    time(imaxPedA2)],     [minSignalCut, maxSignalCut]); lastline('color','g');
line([time(iminPedQ),     time(iminPedQ)],      [minSignalCut, maxSignalCut]); lastline('color','k'); lastline('LineStyle','--');
line([time(iminPedQ2),    time(iminPedQ2)],     [minSignalCut, maxSignalCut]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ),     time(imaxPedQ)],      [minSignalCut, maxSignalCut]); lastline('color','k'); lastline('LineStyle','--');
line([time(imaxPedQ2),    time(imaxPedQ2)],     [minSignalCut, maxSignalCut]); lastline('color','k'); lastline('LineStyle','--');

title('SP wvf (average)'); legend('all','above charge pedestal');
xlabel('time[ns]'); ylabel('A[mV]'); yaxis(minSignalCut, maxSignalCut); xaxis(time(1), time(length(time))); box;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VI. 'SUMMARY'                                                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('                                         ');
disp('               SUMMARY                   ');
disp('                                         ');
disp(['    <nph>         = ', num2str(parQ(1))]);
disp(['    P_0           = ', num2str(exp(-parQ(1)))]);
disp(['    <Q1>          = ', num2str(parQ(4)), ' pC']);
disp(['    sigma_Q1/<Q1> = ', num2str(parQ(3)/parQ(4))]);
if(isAFit)                                                                  %#ok<BDLGI>
    disp(['    <nph> (A)     = ', num2str(parA(1))]);
    disp(['    P_0 (A)       = ', num2str(exp(-parA(1)))]);
    disp(['    <A1>          = ', num2str(parA(4)), ' mV']);
    disp(['    sigma_A1/<A1> = ', num2str(parA(3)/parA(4))]);
end
disp('                                         ');
disp(['    sigma_QPED (single-photon) = ', num2str(parQ(5)),    ' pC']);
disp(['    mean_QPED  (single-photon) = ', num2str(parQ(6)),    ' pC']);
disp(['    sigma_QPED (pure noise)    = ', num2str(parQped(2)), ' pC']);
disp(['    mean_QPED  (pure noise)    = ', num2str(parQped(1)), ' pC']);
if(isAFit) 
disp('                                         ');
disp(['    sigma_APED (single-photon) = ', num2str(parA(5)),    ' mV']);
disp(['    mean_APED  (single-photon) = ', num2str(parA(6)),    ' mV']);
disp(['    sigma_APED (pure noise)    = ', num2str(parAped(2)), ' mV']);
disp(['    mean_APED  (pure noise)    = ', num2str(parAped(1)), ' mV']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% VI. 'SAVE TO FILE'                                                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (isSaveToFile) %#ok<BDLGI>
    fid   = fopen([DIR_OUTPUT, FILE, '_ASCII'], 'w');
    fprintf(fid, '%s\r', 'Nevt   time   V');
    for i=1:Nevts
        if(rem(i, 100)==0), disp('i= ', num2str(i));  end
        for j=1:length(time)
            fprintf(fid, '%d, %3.3f, %3.3f \r', i, time(j), SPsignal(j,i));
        end
    end
end

fclose('all');

mosaic;

% L372, L360, L402