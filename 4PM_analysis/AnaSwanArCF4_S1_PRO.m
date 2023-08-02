%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script for reading data from Tektronix oscilloscope or CAEN DAQ %%%%%
%%%%%                  and writing to ASCII file                      %%%%%
%%%%%                    (DGD 10/Dec/2020)                            %%%%%
%%%%%                  (upgraded 03/Apr/2021)                         %%%%%
%%%%%                   (tuned 12/Jan/2023)                           %%%%%
%%%%%                   (tuned 11/Feb/2023)                           %%%%%
%%%%%                                                                 %%%%%
%%%%% NOTES:                                                          %%%%%
%%%%%                                                                 %%%%%
%%%%% -The script is made for Ar, Xe, Ar/CF4 analysis, relying on PM1 %%%%%
%%%%% and PM4 for trigger and event selection, with PM1/PM4 fit to a  %%%%%
%%%%% to a 2-exponential and PM2/PM3 to a single one. It should be    %%%%%
%%%%% possible to make it more general                                %%%%%
%%%%%                                                                 %%%%%
%%%%% -The isPM3Fit flag seems obsolete, as it is always set to 1.    %%%%%
%%%%%                                                                 %%%%%
%%%%% MAIN DATA:                                                      %%%%%
%%%%% data:    raw waveforms (NPMT, Ntime, Nevt)                      %%%%%
%%%%% data_:   pedestal-subtracted, bin2time, bin2amplitude corrected %%%%%
%%%%% data__:  + cut & time-shifted to the beginning of time-range    %%%%%
%%%%% data___: + calibrated to phe/ns                                 %%%%%
%%%%% data____:+ remove afterpulsing (not implemented yet)            %%%%%
%%%%%                                                                 %%%%%
%%%%% MEAN WVF and STANDARD DEVIATION                                 %%%%%
%%%%% wvf1-4:  pedestal-subtracted, bin2time, bin2amplitude corrected %%%%%
%%%%%          cut, time-shifted, calibrated to phe/ns                %%%%%
%%%%% wvf1-4_: + afterpulsing-subtracted                              %%%%%
%%%%% swvf1-4, swvf1-4_: corresponding spreads                        %%%%%
%%%%%                                                                 %%%%%
%%%%% WORKFLOW:                                                       %%%%%
%%%%%                                                                 %%%%%
%%%%% I.    INITIALIZATION (adjust your desired values here)          %%%%%
%%%%% II.   READOUT (from now on things are automatic -no need to touch) %%
%%%%% III.  START ANALYSIS                                            %%%%%
%%%%%       a) DEFINE WVF WINDOW AND TRIGGER TIME IF NOT AUTOMATIC    %%%%%
%%%%%       b) SUBTRACT PEDESTAL TO WVF                               %%%%%
%%%%%       c) PLOT WVFS FOR COARSE ASSESSMENT                        %%%%%
%%%%% IV.   DIGITIZATION/PULSE SHAPE ANALYSIS                         %%%%%
%%%%% V.    CUTS                                                      %%%%%
%%%%%       a) MAIN                                                   %%%%%
%%%%%       b) STUDY MEAN WVF AFTER CUT                               %%%%%
%%%%%       c) STUDY MEAN WVF AFTER CUT & SHIFT                       %%%%% 
%%%%% VI.   SINGLE-PHOTON ANALYSIS (GENERALLY NOT DONE ->BCKP!)       %%%%%
%%%%% VII.  FIT                                                       %%%%%
%%%%% VIII. CALCULATE W                                               %%%%%
%%%%% IX.   SUMMARY                                                   %%%%% 
%%%%%  X.   SAVE TO FILE                                              %%%%%
%%%%%                                                                 %%%%%
%%%%% TODO:                                                           %%%%%
%%%%%                                                                 %%%%%
%%%%% check pedestal substraction and consider doing wvf by wvf       %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I. INITIALIZATION (adjust your desired values here)                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

NevtsAna = 0;                                                               %Maximum number of events to analyze. Set to zero if not to be used.

% READ
DataType = 'CAEN';   % 'CAEN' // 'OSCI'                                     %Specify data type to be read (CAEN card and Tektronix osci, so far)
Nsave    = 20000;                                                           %Number of events from card/scope to be saved in workspace (workspace limits in Matlab 8)
nCh      = 4;                                                               %Number of wvf channels to be read
FILE     = 'wave';                                                          %Name of file to be read and saved

% ANA
pedEntries   = 40;                                                          %Number of entries to be used for pedestal estimate
isSmearedA   = 1;                                                           %Amplitude smearing within a bin
isSmearedT   = 1;                                                           %Time smearing within a bin
NPMTs        = 4;                                                           %Number of PMT wvf recorded
isPMforNorm  = 1;                                                           %PM for relative normalization (e.g., to estimate position, default is 1)
isPMforShift = 1;                                                           %PM for relative shift (default is 1)
Rin          = 50;                                                          %channel input impedance
sigmaVth     = 3;                                                           %Number of sigmas above threshold for pulse shape analysis (in pple do not touch!)
                                                                            %NOTE: use common threshold for all wvfs per PM (may be done individually in the future)
% TRIGGER
isTriggerSetAutomatically = 1;                                              %Find trigger region from average wvf in reference channel (default is PM1)
AmplitudeWindow = [500, 1500];                                              %If automatic does not work, chose manual

% CUTS
isA1cut  = 1*0; isA4cut = 1*0;
isXcut   = 1*0; isYcut = 1*0;                                               %Apply position cut in X (PM2, PM3) or Y (PM1, PM4)
isT20cut = 1;                                                               %Apply cut in time approximately where pulse 'starts' (defined at 20% of signal maximum)

% FIT
Tmax1 = 180; Tmax2 = 50; Tmax3 = 50; Tmax4 = 180;                           %Time range for fit (excludes afterpulsing and noise)
isPM3Fit = 1;                                                               %Temporary, allows to deal with some cases where the PM3 fit does not work well.
                                                                            %PM1 and PM4 are fit to double exponentials and PM2 and PM3 to single. NOTE: it might be made more general.
% SINGLE-PHOTON CALIBRATION
isSPanalysis = 1*0;                                                         %Do Single-Photon analysis with the signal tails (requires additional tuning)
DeltaT_PM1 = 20;  DeltaT_PM2 = 20;  DeltaT_PM3 = 20;  DeltaT_PM4 = 20;      %Time window
Toff_PM1   = 160; Toff_PM2   = 15;  Toff_PM3   = 20;  Toff_PM4   = 160;     %Delay relative to trigger 
sumA_to_ph = 5;                                                             %Convert to photons in approximate way (temporary, needs fix)                 

% AUXILLIARY PARAMETERS (for representation, etc)
Arange     = -10:1:100;                                                     %Range for Amplitude representation [mV]
Qrange     = -1:0.1:10;                                                     %Range for Charge representation [nC]

% INITIALIZATION
AnaSwanArCF4_S1ini;                                                         %Load input file and over-write some of the above parameters if needed, case by case.

% SAVE
isSaveSummaryToFile = 1;                                                    %Flag to indicate whether summary saved to file or not
isSaveWVFsToFile    = 1;                                                    %Flag to indicate whether wvfs    saved to file or not
NwvfsWritten        = 201;                                                  %Number of waveforms to be written to file (if 0 write all)
DIR_OUTPUT          = DIR;                                                  %Only if something is saved

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% II. READOUT (from now on things are automatic -no need to touch)    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(DataType,'OSCI')), ReadDataOsci; end
if(strcmp(DataType,'CAEN')), ReadDataCAEN_PRO; end

if(NevtsAna)
    Nevts = NevtsAna;
    data = data(:, :, 1:Nevts);
end

%%%%%%%%%%%%%%%%%%%%
%III. START ANALYSIS
%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%III.a) DEFINE WVF WINDOW AND TRIGGER TIME IF NOT AUTOMATIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AmpMeanWVF = zeros(1,NPMTs); TimeMeanWVF = zeros(1,NPMTs); tMinWVF = zeros(1,NPMTs); tMaxWVF = zeros(1,NPMTs);
 
if(isTriggerSetAutomatically)
    disp('note: trigger window will be set automatically!');
    for k=1:NPMTs
        MeanWVF             = mean(squeeze(data(k,:,:))'); %#ok<UDIM>
        PedMeanWVF          = MeanWVF(1);
        MeanWVF             = MeanWVF - PedMeanWVF;
        AmpMeanWVF(k)       = max(MeanWVF);
        TimeMeanWVF(k)      = time(MeanWVF==AmpMeanWVF(k));
        
        [~, ~, tMinWVF(k), tMaxWVF(k), ~, ~, ~, ~, ~, ~, ~, ~] = PSA_Nausicaa0(MeanWVF, time, AmpMeanWVF(k)/10, TimeMeanWVF(k));
    end
    trigTime      = mean(TimeMeanWVF);
    trigWindow(1) = min(tMinWVF);
    trigWindow(2) = min(tMaxWVF);
else
    trigTime      = 0.5*(AmplitudeWindow(2) + AmplitudeWindow(1));
    trigWindow(1) = AmplitudeWindow(1);
    trigWindow(2) = AmplitudeWindow(2);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%III.b) SUBTRACT PEDESTAL TO WVF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pedMEAN = zeros(NPMTs, Nevts); pedSTD = zeros(NPMTs, Nevts);      

data_  = zeros(size(data));                                                 %'_' == pedestal-subtracted data

for k=1:NPMTs
    for i=1:Nevts
        dataTMP         = data(k,:,i);
        data_ped(k,:,i) = dataTMP(time<pedEntries*TBin);                    %#ok<SAGROW>
        
        pedMEAN(k,i)    = mean(data_ped(k,:,i));                            %used for pedestal subtraction
        pedSTD(k,i)     = std(data_ped(k,:,i));                             %used for threshold determination in PSA
        
        data_(k,:,i)    = data(k,:,i) - pedMEAN(k,i);                       %pedestal subtraction
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%III.c) PLOT WVFS FOR COARSE ASSESSMENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
for k=1:NPMTs
    dataMean(k,:) = mean(squeeze(data_(k,:,:))');                           %#ok<UDIM,SAGROW>
    subplot(NPMTs,1,k);hold on;
    plot(time, data_(k,:,1)); plot(time, data_(k,:,2), 'r'); plot(time, data_(k,:,3), 'g'); plot(time, data_(k,:,4), 'c'); plot(time, data_(k,:,5), 'k');
    maxYForWVF = mean(max(data_(k,:,:)))*1.5; 
    minYForWVF = -2;                      %Ymin for representation
    plot(time, dataMean(k,:),'k'); lastline('LineWidth',3);
    line([trigWindow(1), trigWindow(1)], [-500,500]); lastline('color','k');  lastline('LineStyle','--');
    line([trigWindow(2), trigWindow(2)], [-500 500]); lastline('color','k');  lastline('LineStyle','--');    
    line([trigTime,      trigTime],      [-500 500]); lastline('color','r'); 
    if(k==1)
        title('PMT 1 wvfs');
        legend('wvf 1', 'wvf 2', 'wvf 3', 'wvf 4', 'wvf 5', 'Mean all evt-wvf');
    end
    if(k==2),title('PMT 2 wvfs');end
    if(k==3),title('PMT 3 wvfs');end
    if(k==4),title('PMT 4 wvfs');end
    xaxis(min(time),max(time)); yaxis(minYForWVF, maxYForWVF); 
    xlabel('time[ns]'); ylabel('A[mV]'); box on;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% IV. 'DIGITIZATION/PULSE SHAPE ANALYSIS': %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Amp       = zeros(NPMTs, Nevts); charge    = zeros(NPMTs, Nevts); meanT = zeros(NPMTs, Nevts); 
trise_old = zeros(NPMTs, Nevts); tMin      = zeros(NPMTs, Nevts); trise = zeros(NPMTs, Nevts);
rmsOUT    = zeros(NPMTs, Nevts); stdT      = zeros(NPMTs, Nevts); Vth   = zeros(NPMTs, 1);

for k=1:NPMTs
    Vth(k) = sigmaVth*mean(pedSTD(k,:));
    
    for i=1:Nevts
        F = data_(k,:,i);  t = time;
        [charge(k,i), Amp(k,i), tMin(k,i), ~, ~, meanT(k,i), trise_old(k,i), ~, rmsOUT(k,i), trise(k,i), ~, ~, ~, ~] = PSA_Nausicaa0(F, t, Vth(k), trigTime);
        
        if(isSmearedA), Amp(k,i)    = Amp(k,i) - ABin/2 + ABin*rand; end
        
        charge(k,i)    = charge(k,i) * TBin / Rin;                   %[mV x ns = pC]
    end
    disp(['PM ', num2str(k), ' digitized']);
end

%Define PM and event -based magnitudes
A1      = Amp(1,:);       A2      = Amp(2,:);     A3      = Amp(3,:);     A4      = Amp(4,:);
Q1      = charge(1,:);    Q2      = charge(2,:);  Q3      = charge(3,:);  Q4      = charge(4,:);
Q1norm  = Q1*mean(isPMforNorm)/mean(Q1);  Q2norm  = Q2*mean(isPMforNorm)/mean(Q2);
Q3norm  = Q3*mean(isPMforNorm)/mean(Q3);  Q4norm  = Q4*mean(isPMforNorm)/mean(Q4);  

T1      = meanT(1,:);     T2      = meanT(2,:);   T3      = meanT(3,:);   T4      = meanT(4,:);
stdP1   = pedSTD(1,:);    stdP2   = pedSTD(2,:);  stdP3   = pedSTD(3,:);  stdP4   = pedSTD(4,:);
Base1   = rmsOUT(1,:);    Base2   = rmsOUT(2,:);  Base3   = rmsOUT(3,:);  Base4   = rmsOUT(4,:);

T20_1   = trise_old(1,:) + tMin(1,:) - trise(1,:);                          %Trick to get from the PSA-function the risetime at 20% level.
T20_2   = trise_old(2,:) + tMin(2,:) - trise(2,:);
T20_3   = trise_old(3,:) + tMin(3,:) - trise(3,:);
T20_4   = trise_old(4,:) + tMin(4,:) - trise(4,:);

Yevt  = zeros(1,Nevts); Xevt  = zeros(1,Nevts); Revt  = zeros(1,Nevts);

for i=1:Nevts   
    Yevt(i) = (Q1norm(i)/(Q1norm(i) + Q4norm(i))-0.5) * 2;
    Xevt(i) = (Q3norm(i)/(Q3norm(i) + Q2norm(i))-0.5) * 2;
    Revt(i) = sqrt(Xevt(i)*Xevt(i) + Yevt(i)*Yevt(i));
end

%PLOT X,Y,R distributions
figure;
subplot(4,1,1); hold on;
for k=1:NPMTs
    if(k==1), hist1D(Q1norm, 0:0.01:5); lastline('color', 'k');end
    if(k==2), hist1D(Q2norm, 0:0.01:5); lastline('color', 'r');end
    if(k==3), hist1D(Q3norm, 0:0.01:5); lastline('color', 'g');end
    if(k==4), hist1D(Q4norm, 0:0.01:5); lastline('color', 'b');end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('Q[pC]'); ylabel('entries'); title ('Charged (normalized to ref PM)'); box on;
subplot(4,1,2);
plot(Xevt,Yevt,'.'); lastline('MarkerSize', 1);
xlabel('x_{evt} [mm]'); ylabel('y_{evt} [mm]'); title('position before cuts');
xaxis(-2, 2); yaxis(-2, 2);
subplot(4,1,3)
hist1D(Xevt,-2:0.01:2);
xlabel('x_{evt} (3-2) [mm]');
subplot(4,1,4)
hist1D(Yevt, -2:0.01:2);
xlabel('y_{evt} (4-1) [mm]');

%PLOT Q-Q correlations
figure; 
subplot(2,2,1);
plot(Q1, Q2, '.'); title('PM1 vs PM2');
xlabel('Q_{1}[pC]'); ylabel('Q_{2}[pC]');
subplot(2,2,2);
plot(Q1, Q3, '.'); title('PM1 vs PM3');
xlabel('Q_{1}[pC]'); ylabel('Q_{3}[pC]');
subplot(2,2,3);
plot(Q1, Q4, '.'); title('PM1 vs PM4');
xlabel('Q_{1}[pC]'); ylabel('Q_{4}[pC]');
subplot(2,2,4);
plot(Q2, Q4, '.'); title('PM2 vs PM4');
xlabel('Q_{2}[pC]'); ylabel('Q_{4}[pC]'); 

%PLOT PM variables
figure;
subplot(2,1,1);
hold on;

for k=1:NPMTs
    if(k==1), hist1D(A1, Arange); lastline('color', 'k');end
    if(k==2), hist1D(A2, Arange); lastline('color', 'r');end
    if(k==3), hist1D(A3, Arange); lastline('color', 'g');end
    if(k==4), hist1D(A4, Arange); lastline('color', 'b');end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xaxis(min(Arange), max(Arange));
xlabel('A[mV]'); ylabel('entries'); title ('Amplitude'); box on;

subplot(2,1,2);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(Q1, Qrange); lastline('color', 'k');end
    if(k==2), hist1D(Q2, Qrange); lastline('color', 'r');end
    if(k==3), hist1D(Q3, Qrange); lastline('color', 'g');end
    if(k==4), hist1D(Q4, Qrange); lastline('color', 'b');end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xaxis(min(Qrange), max(Qrange));
xlabel('Q[pC]'); ylabel('entries'); title ('Charge'); box on;

figure;
subplot(1,2,1);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(stdP1, -0.1:0.05:5); lastline('color', 'k'); end
    if(k==2), hist1D(stdP2, -0.1:0.05:5); lastline('color', 'r'); end
    if(k==3), hist1D(stdP3, -0.1:0.05:5); lastline('color', 'g'); end
    if(k==4), hist1D(stdP4, -0.1:0.05:5); lastline('color', 'b'); end
end
for k=1:NPMTs
    if(k==1), line([Vth(k), Vth(k)], [0,2000]); lastline('color', 'k'); end
    if(k==2), line([Vth(k), Vth(k)], [0,2000]); lastline('color', 'r'); end
    if(k==3), line([Vth(k), Vth(k)], [0,2000]); lastline('color', 'g'); end
    if(k==4), line([Vth(k), Vth(k)], [0,2000]); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xaxis(0,5);
xlabel('A[mV]'); ylabel('entries'); title ('pedestal std (pre wvf)'); box on;

subplot(1,2,2);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(Base1, -0.1:0.05:5); lastline('color', 'k'); end
    if(k==2), hist1D(Base2, -0.1:0.05:5); lastline('color', 'r'); end
    if(k==3), hist1D(Base3, -0.1:0.05:5); lastline('color', 'g'); end
    if(k==4), hist1D(Base4, -0.1:0.05:5); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xaxis(-0.1,5);
xlabel('base fluct [mV]'); ylabel('entries'); title ('pedestal std (pre+post wvf)'); box on;

figure;
hold on;
for k=1:NPMTs
    if(k==1), hist1D(T1, trigWindow(1):trigWindow(2)); lastline('color', 'k'); end
    if(k==2), hist1D(T2, trigWindow(1):trigWindow(2)); lastline('color', 'r'); end
    if(k==3), hist1D(T3, trigWindow(1):trigWindow(2)); lastline('color', 'g'); end
    if(k==4), hist1D(T4, trigWindow(1):trigWindow(2)); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('t[ns]'); ylabel('entries'); title ('mean time'); box on;

figure;
hold on;
for k=1:NPMTs
    if(k==1), hist1D(T20_1, trigWindow(1):trigWindow(2)); lastline('color', 'k'); end
    if(k==2), hist1D(T20_2, trigWindow(1):trigWindow(2)); lastline('color', 'r'); end
    if(k==3), hist1D(T20_3, trigWindow(1):trigWindow(2)); lastline('color', 'g'); end
    if(k==4), hist1D(T20_4, trigWindow(1):trigWindow(2)); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('t[ns]'); ylabel('entries'); title ('time at 20%'); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V. 'CUTS': %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  V.a): MAIN  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Icut   = A1>0;                                                              %No cut. All subsequent cuts are added to this
if(isA1cut),  A1_LOW = 3*mean(stdP1); Icut = Icut & A1>A1_LOW;    end %#ok<BDLGI>
if(isA4cut),  A4_LOW = 3*mean(stdP4); Icut = Icut & A4>A4_LOW;    end %#ok<BDLGI>

if(isT20cut)
    T20_l_tmp = T20_1(find((T20_1>=trigWindow(1)) & (T20_1<=trigWindow(2))));
    [N,X]     = hist1D(T20_l_tmp,trigWindow(1):TBin:trigWindow(2));
    Xmax      = mean(X(find(N==max(N))));
    T20_UP    = Xmax + 3*std(T20_l_tmp);
    T20_LOW   = Xmax - 3*std(T20_l_tmp);
    Icut      = Icut & T20_1<T20_UP & T20_1>T20_LOW;
end
if(isXcut) %#ok<BDLGI>
    Icut = Icut & Xevt>-2 & Xevt<2;                                         %Removes outliers
    rmsX = std(Xevt(Icut)); meanX = mean(Xevt(Icut));
    Xmin = meanX - 1.5*rmsX; Xmax = meanX + 1.5*rmsX;
    Icut = Icut & (Xevt>Xmin) & (Xevt<Xmax); 
end
if(isYcut) %#ok<BDLGI>
    Icut = Icut & Yevt>-2 & Yevt<2;                                         %Removes outliers   
    rmsY = std(Yevt(Icut)); meanY = mean(Yevt(Icut));
    Ymin = meanY - 1.5*rmsY; Ymax = meanY + 1.5*rmsY;
    Icut = Icut & (Yevt>Ymin) & (Yevt<Ymax);  
end

%PLOT X,Y,R distributions (after cuts)
figure;
subplot(4,1,1)
plot(Xevt(Icut),Yevt(Icut),'.'); lastline('MarkerSize', 1);
xlabel('x_{evt} [mm]'); ylabel('y_{evt}[mm]'); title('position after cuts [cut]');
xaxis(-2, 2); yaxis(-2, 2);
subplot(4,1,2)
hist1D(Revt(Icut), 0:0.01:2);
xlabel('R_{evt} [mm]');
subplot(4,1,3)
hist1D(Xevt(Icut), -2:0.01:2);
if(isXcut) %#ok<BDLGI>
    line([Xmin, Xmin], [0, max(hist1D(Xevt(Icut), -2:0.01:2))]); lastline('color','r');  lastline('LineStyle','--');
    line([Xmax, Xmax], [0, max(hist1D(Xevt(Icut), -2:0.01:2))]); lastline('color','r');  lastline('LineStyle','--');
end
xlabel('X_{evt} [mm]');
subplot(4,1,4)
hist1D(Yevt(Icut), -2:0.01:2);
if(isYcut) %#ok<BDLGI>
    line([Ymin, Ymin], [0, max(hist1D(Yevt(Icut), -2:0.01:2))]); lastline('color','r');  lastline('LineStyle','--');
    line([Ymax, Ymax], [0, max(hist1D(Yevt(Icut), -2:0.01:2))]); lastline('color','r');  lastline('LineStyle','--');
end
xlabel('Y_{evt} [mm]');

%PLOT Q-Q correlations (after cuts)
figure; 
subplot(2,2,1);
plot(Q1(Icut), Q2(Icut), '.'); title('PM1 vs PM2');
xlabel('Q_{1}[pC]'); ylabel('Q_{2}[pC]');
subplot(2,2,2);
plot(Q1(Icut), Q3(Icut), '.'); title('PM1 vs PM3');
xlabel('Q_{1}[pC]'); ylabel('Q_{3}[pC]');
subplot(2,2,3);
plot(Q1(Icut), Q4(Icut), '.'); title('PM1 vs PM4');
xlabel('Q_{1}[pC]'); ylabel('Q_{4}[pC]');
subplot(2,2,4);
plot(Q2(Icut), Q4(Icut), '.'); title('PM2 vs PM4');
xlabel('Q_{2}[pC]'); ylabel('Q_{4}[pC]'); 

%PLOT PM variables (after cuts)
figure;
subplot(2,1,1);
hold on;

for k=1:NPMTs
    if(k==1), hist1D(A1(Icut), Arange); lastline('color', 'k');end
    if(k==2), hist1D(A2(Icut), Arange); lastline('color', 'r');end
    if(k==3), hist1D(A3(Icut), Arange); lastline('color', 'g');end
    if(k==4), hist1D(A4(Icut), Arange); lastline('color', 'b');end
end
if(isA1cut) %#ok<BDLGI>
    line([A1_LOW, A1_LOW], [0, max(hist1D(A1(Icut), Arange))]); lastline('color','k');  lastline('LineStyle','--');
end
if(isA4cut) %#ok<BDLGI>
    line([A4_LOW, A4_LOW], [0, max(hist1D(A4(Icut), Arange))]); lastline('color','b');  lastline('LineStyle','--');
end
xaxis(min(Arange), max(Arange));
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('A[mV]'); ylabel('entries'); title ('Amplitude [cut]'); box on;

subplot(2,1,2);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(Q1(Icut), Qrange); lastline('color', 'k');end
    if(k==2), hist1D(Q2(Icut), Qrange); lastline('color', 'r');end
    if(k==3), hist1D(Q3(Icut), Qrange); lastline('color', 'g');end
    if(k==4), hist1D(Q4(Icut), Qrange); lastline('color', 'b');end
end
xaxis(min(Qrange), max(Qrange));
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('Q[pC]'); ylabel('entries'); title ('Charge [cut]'); box on;

figure;
subplot(1,2,1);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(stdP1(Icut), -0.1:0.05:5); lastline('color', 'k'); end
    if(k==2), hist1D(stdP2(Icut), -0.1:0.05:5); lastline('color', 'r'); end
    if(k==3), hist1D(stdP3(Icut), -0.1:0.05:5); lastline('color', 'g'); end
    if(k==4), hist1D(stdP4(Icut), -0.1:0.05:5); lastline('color', 'b'); end
end
N1 = max(hist1D(stdP1(Icut), -0.1:0.05:5)); N2 = max(hist1D(stdP2(Icut), -0.1:0.05:5));
N3 = max(hist1D(stdP3(Icut), -0.1:0.05:5)); N4 = max(hist1D(stdP4(Icut), -0.1:0.05:5));

for k=1:NPMTs
    if(k==1), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'k'); lastline('LineStyle', '--'); end
    if(k==2), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'r'); lastline('LineStyle', '--'); end
    if(k==3), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'g'); lastline('LineStyle', '--'); end
    if(k==4), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'b'); lastline('LineStyle', '--'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4', 'Vth1', 'Vth2', 'Vth3', 'Vth4');
xaxis(0,5);
xlabel('A[mV]'); ylabel('entries'); title ('pedestal std (pre wvf) [cut]'); box on;

subplot(1,2,2);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(Base1(Icut), -0.1:0.05:5); lastline('color', 'k'); end
    if(k==2), hist1D(Base2(Icut), -0.1:0.05:5); lastline('color', 'r'); end
    if(k==3), hist1D(Base3(Icut), -0.1:0.05:5); lastline('color', 'g'); end
    if(k==4), hist1D(Base4(Icut), -0.1:0.05:5); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xaxis(-0.1,5);
xlabel('base fluct [mV]'); ylabel('entries'); title ('pedestal std (pre+post wvf) [cut]'); box on;

figure;
hold on;
for k=1:NPMTs
    if(k==1), hist1D(T1(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'k'); end
    if(k==2), hist1D(T2(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'r'); end
    if(k==3), hist1D(T3(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'g'); end
    if(k==4), hist1D(T4(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'b'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('t[ns]'); ylabel('entries'); title ('mean time [cut]'); box on;

figure;
hold on;
for k=1:NPMTs
    if(k==1), hist1D(T20_1(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'k'); end
    if(k==2), hist1D(T20_2(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'r'); end
    if(k==3), hist1D(T20_3(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'g'); end
    if(k==4), hist1D(T20_4(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'b'); end
end
if(isT20cut)
    line([T20_LOW, T20_LOW], [0, max(hist1D(T20_1(Icut), trigWindow(1):trigWindow(2)))]); lastline('color','k');  lastline('LineStyle','--');
    line([T20_UP,  T20_UP],  [0, max(hist1D(T20_1(Icut), trigWindow(1):trigWindow(2)))]); lastline('color','k');  lastline('LineStyle','--');
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
xlabel('t[ns]'); ylabel('entries'); title ('time at 20% [cut]'); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  V.b) STUDY MEAN WVF AFTER CUT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PMT_MWVF_cut = zeros(NPMTs, dataSize(2));
PMT_MWVF     = zeros(NPMTs, dataSize(2));
for k=1:NPMTs
    for i=1:Nevts
        if(Icut(i)),             PMT_MWVF_cut(k,:) = PMT_MWVF_cut(k,:) + data_(k,:,i); end
        PMT_MWVF(k,:) = PMT_MWVF(k,:) + data_(k,:,i);
    end
end
PMT_MWVF_cut = PMT_MWVF_cut/length(find(Icut));
PMT_MWVF     = PMT_MWVF/length(Icut);

figure; 
subplot(2,2,1); hold on;
plot(time, PMT_MWVF(1,:));     lastline('color','k');   lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(1,:)); lastline('color','k');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs no-cut)'); box on;

subplot(2,2,2); hold on;
plot(time, PMT_MWVF(2,:));     lastline('color','b');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(2,:)); lastline('color','b');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs no-cut)'); box on;

subplot(2,2,3); hold on;
plot(time, PMT_MWVF(3,:));     lastline('color','r');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(3,:)); lastline('color','r');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs no-cut)'); box on;

subplot(2,2,4); hold on;
plot(time, PMT_MWVF(4,:));     lastline('color','g');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(4,:)); lastline('color','g');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs no-cut)'); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  V.c): STUDY MEAN WVF AFTER CUT & SHIFT  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PMT_MWVF_cut_shifted = zeros(NPMTs, dataSize(2));                           %Mean
PMT_SWVF_cut_shifted = zeros(NPMTs, dataSize(2));                           %Error of the mean
data__ = zeros(size(data));
data__ = data__(:, :, 1:length(A1(Icut)));

%(1) Shift ALL wvfs
Ishift   = 4;                                                               %Bin number where to realign waveform
if(Ishift>pedEntries), Ishift = pedEntries; end                             %Do not make larger than pedestal (NOTE: really needed??)
for k=1:NPMTs
    ievt = 0;
    for i=1:Nevts
        if(Icut(i))
            ievt = ievt + 1;
            tf = time(Ishift); ti = T20_1(i);
            wvf_i = squeeze(data_(k,:,i));
            wvf_f = TimeAlign(time, wvf_i, ti, tf);
            PMT_MWVF_cut_shifted(k,:) = PMT_MWVF_cut_shifted(k,:) + wvf_f;
            PMT_SWVF_cut_shifted(k,:) = PMT_SWVF_cut_shifted(k,:) + wvf_f.^2;
            data__(k, :, ievt) = wvf_f;
        end
    end
end

for k=1:NPMTs
    PMT_MWVF_cut_shifted(k,:) = PMT_MWVF_cut_shifted(k,:)/length(find(Icut));
    PMT_SWVF_cut_shifted(k,:) = sqrt(PMT_SWVF_cut_shifted(k,:)/length(find(Icut)) - PMT_MWVF_cut_shifted(k,:).*PMT_MWVF_cut_shifted(k,:));
    PMT_SWVF_cut_shifted(k,:) = PMT_SWVF_cut_shifted(k,:)/sqrt(length(find(Icut)));
end

%(2) Shift average waveforms to PM1 maximum (approximate). NOTE: PM1 by default, does not matter except numerically perhaps
for k=1:NPMTs    
    Iref = find(PMT_MWVF_cut(isPMforShift,:)==max(PMT_MWVF_cut(isPMforShift,:)));
    I = find(PMT_MWVF_cut(k,:)==max(PMT_MWVF_cut(k,:)));
    tf = time(Iref); ti = time(I);
    PMT_MWVF_cut(k, :) = TimeAlign(time, PMT_MWVF_cut(k, :), ti, tf);
    
    Iref = find(PMT_MWVF_cut_shifted(isPMforShift,:)==max(PMT_MWVF_cut_shifted(isPMforShift,:)));
    I = find(PMT_MWVF_cut_shifted(k,:)==max(PMT_MWVF_cut_shifted(k,:)));
    tf = time(Iref); ti = time(I);
    PMT_MWVF_cut_shifted(k, :) = TimeAlign(time, PMT_MWVF_cut_shifted(k, :), ti, tf);
    wvf = TimeAlign(time, PMT_SWVF_cut_shifted(k, :), ti, tf);
    wvf(find(wvf==0)) = max(PMT_SWVF_cut_shifted(k, :));                    %Reassign error to zeros for later fit
    PMT_SWVF_cut_shifted(k, :) = wvf;    
end

%(3) Shift average waveforms to leave only Ishift bins from 0 to max
for k=1:NPMTs    
    I = find(PMT_MWVF_cut(k,:)==max(PMT_MWVF_cut(k,:)));
    tf = time(Ishift); ti = time(I);
    PMT_MWVF_cut(k, :) = TimeAlign(time, PMT_MWVF_cut(k, :), ti, tf);
    
    I = find(PMT_MWVF_cut_shifted(k,:)==max(PMT_MWVF_cut_shifted(k,:)));
    tf = time(Ishift); ti = time(I);
    PMT_MWVF_cut_shifted(k, :) = TimeAlign(time, PMT_MWVF_cut_shifted(k, :), ti, tf);
    wvf = TimeAlign(time, PMT_SWVF_cut_shifted(k, :), ti, tf);
    wvf(find(wvf==0)) = max(PMT_SWVF_cut_shifted(k, :));                   
    PMT_SWVF_cut_shifted(k, :) = wvf;    
end

%(4) Re-plot average waveforms (both cut & shifted) by selecting regions 
% with S/N larger than 10
% NOTE: it is generally not used and could be moved to function
% (StoNcut_PM is usually defined as 0 in the initialization file)
thres_MWVF1 = StoNcut_PM1*std(PMT_MWVF_cut(1,1:pedEntries));
thres_MWVF2 = StoNcut_PM2*std(PMT_MWVF_cut(2,1:pedEntries));
thres_MWVF3 = StoNcut_PM3*std(PMT_MWVF_cut(3,1:pedEntries));
thres_MWVF4 = StoNcut_PM4*std(PMT_MWVF_cut(4,1:pedEntries));

for k=1:NPMTs
    if(k==1 && isPM1_StoNcut)
        PMT_MWVF_cut1 = PMT_MWVF_cut(k, :);
        maxwvf = max(PMT_MWVF_cut1);
        I = find((PMT_MWVF_cut1<thres_MWVF1 & time<mean(time(find(PMT_MWVF_cut1==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut1<thres_MWVF1 & time>mean(time(find(PMT_MWVF_cut1==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut1 = PMT_MWVF_cut1(imin:imax);
        tim_MWVF_cut1 = time(imin:imax);
        
        PMT_MWVF_cut_shifted1 = PMT_MWVF_cut_shifted(k, :);
        maxwvf = max(PMT_MWVF_cut_shifted1);
        I = find((PMT_MWVF_cut_shifted1<thres_MWVF1 & time<mean(time(find(PMT_MWVF_cut_shifted1==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut_shifted1<thres_MWVF1 & time>mean(time(find(PMT_MWVF_cut_shifted1==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut_shifted1 = PMT_MWVF_cut_shifted1(imin:imax);
        tim_MWVF_cut_shifted1 = time(imin:imax);

        PMT_SWVF_cut_shifted1 = PMT_SWVF_cut_shifted(k, imin:imax);
    elseif(k==1 && ~isPM1_StoNcut)
        PMT_MWVF_cut1 = PMT_MWVF_cut(k, :);        
        tim_MWVF_cut1 = time;        
        PMT_MWVF_cut_shifted1 = PMT_MWVF_cut_shifted(k, :);        
        PMT_SWVF_cut_shifted1 = PMT_SWVF_cut_shifted(k, :);
        tim_MWVF_cut_shifted1 = time;
    end
    if(k==2 && isPM2_StoNcut)        
        PMT_MWVF_cut2 = PMT_MWVF_cut(k, :);
        maxwvf = max(PMT_MWVF_cut2);
        I = find((PMT_MWVF_cut2<thres_MWVF2 & time<mean(time(find(PMT_MWVF_cut2==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut2<thres_MWVF2 & time>mean(time(find(PMT_MWVF_cut2==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut2 = PMT_MWVF_cut2(imin:imax);
        tim_MWVF_cut2 = time(imin:imax);
        
        PMT_MWVF_cut_shifted2 = PMT_MWVF_cut_shifted(k, :);
        maxwvf = max(PMT_MWVF_cut_shifted2);
        I = find((PMT_MWVF_cut_shifted2<thres_MWVF2 & time<mean(time(find(PMT_MWVF_cut_shifted2==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut_shifted2<thres_MWVF2 & time>mean(time(find(PMT_MWVF_cut_shifted2==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut_shifted2 = PMT_MWVF_cut_shifted2(imin:imax);        
        PMT_SWVF_cut_shifted2 = PMT_SWVF_cut_shifted(k, imin:imax);
        tim_MWVF_cut_shifted2 = time(imin:imax);
    elseif(k==2 && ~isPM2_StoNcut)
        PMT_MWVF_cut2 = PMT_MWVF_cut(k, :);        
        tim_MWVF_cut2 = time;        
        PMT_MWVF_cut_shifted2 = PMT_MWVF_cut_shifted(k, :);   
        PMT_SWVF_cut_shifted2 = PMT_SWVF_cut_shifted(k, :);
        tim_MWVF_cut_shifted2 = time;
    end
    if(k==3 && isPM3_StoNcut)
        PMT_MWVF_cut3 = PMT_MWVF_cut(k, :);
        maxwvf = max(PMT_MWVF_cut3);
        I = find((PMT_MWVF_cut3<thres_MWVF3 & time<mean(time(find(PMT_MWVF_cut3==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut3<thres_MWVF3 & time>mean(time(find(PMT_MWVF_cut3==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut3 = PMT_MWVF_cut3(imin:imax);
        tim_MWVF_cut3 = time(imin:imax);
        
        PMT_MWVF_cut_shifted3 = PMT_MWVF_cut_shifted(k, :);
        maxwvf = max(PMT_MWVF_cut_shifted3);
        I = find((PMT_MWVF_cut_shifted3<thres_MWVF3 & time<mean(time(find(PMT_MWVF_cut_shifted3==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut_shifted3<thres_MWVF3 & time>mean(time(find(PMT_MWVF_cut_shifted3==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut_shifted3 = PMT_MWVF_cut_shifted3(imin:imax);
        tim_MWVF_cut_shifted3 = time(imin:imax);

        PMT_SWVF_cut_shifted3 = PMT_SWVF_cut_shifted(k, imin:imax);
    elseif(k==3 && ~isPM3_StoNcut)
        PMT_MWVF_cut3 = PMT_MWVF_cut(k, :);        
        tim_MWVF_cut3 = time;        
        PMT_MWVF_cut_shifted3 = PMT_MWVF_cut_shifted(k, :);        
        PMT_SWVF_cut_shifted3 = PMT_SWVF_cut_shifted(k, :);
        tim_MWVF_cut_shifted3 = time;
    end
    if(k==4 && isPM4_StoNcut)
        PMT_MWVF_cut4 = PMT_MWVF_cut(k, :);
        maxwvf = max(PMT_MWVF_cut4);
        I = find((PMT_MWVF_cut4<thres_MWVF4 & time<mean(time(find(PMT_MWVF_cut4==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut4<thres_MWVF4 & time>mean(time(find(PMT_MWVF_cut4==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut4 = PMT_MWVF_cut4(imin:imax);
        tim_MWVF_cut4 = time(imin:imax);
        
        PMT_MWVF_cut_shifted4 = PMT_MWVF_cut_shifted(k, :);
        maxwvf = max(PMT_MWVF_cut_shifted4);
        I = find((PMT_MWVF_cut_shifted4<thres_MWVF4 & time<mean(time(find(PMT_MWVF_cut_shifted4==maxwvf)))));
        if(~isempty(I)), imin = max(I)+1; 
        else,            imin = 1; 
        end
        I = find((PMT_MWVF_cut_shifted4<thres_MWVF4 & time>mean(time(find(PMT_MWVF_cut_shifted4==maxwvf)))));
        if(~isempty(I)), imax = min(I)-1; 
        else,            imax = length(time); 
        end
        PMT_MWVF_cut_shifted4 = PMT_MWVF_cut_shifted4(imin:imax);
        tim_MWVF_cut_shifted4 = time(imin:imax);

        PMT_SWVF_cut_shifted4 = PMT_SWVF_cut_shifted(k, imin:imax);
    elseif(k==4 && ~isPM4_StoNcut)
        PMT_MWVF_cut4 = PMT_MWVF_cut(k, :);        
        tim_MWVF_cut4 = time;        
        PMT_MWVF_cut_shifted4 = PMT_MWVF_cut_shifted(k, :);        
        PMT_SWVF_cut_shifted4 = PMT_SWVF_cut_shifted(k, :);
        tim_MWVF_cut_shifted4 = time;
    end
end

figure; 
subplot(2,2,1); hold on;
plot(tim_MWVF_cut_shifted1, PMT_MWVF_cut_shifted1); lastline('color','k');
legend('time-aligned evt by evt');
if(isPM1_StoNcut)
    plot(tim_MWVF_cut1,         PMT_MWVF_cut1);         lastline('color','k');  lastline('LineStyle','--');
    legend('time-aligned evt by evt', 'no time-aligned'); logy; box;
end
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (shifted & cut)'); box on;

subplot(2,2,2); hold on;
plot(tim_MWVF_cut_shifted2, PMT_MWVF_cut_shifted2); lastline('color','b');
legend('time-aligned evt by evt');
if(isPM2_StoNcut)
    plot(tim_MWVF_cut2,         PMT_MWVF_cut2);         lastline('color','b');  lastline('LineStyle','--');
    legend('time-aligned evt by evt', 'no time-aligned'); logy; box;
end
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (shifted & cut)'); box on;

subplot(2,2,3); hold on;
plot(tim_MWVF_cut_shifted3, PMT_MWVF_cut_shifted3); lastline('color','r');
legend('time-aligned evt by evt');
if(isPM3_StoNcut)
    plot(tim_MWVF_cut3,         PMT_MWVF_cut3);         lastline('color','r');  lastline('LineStyle','--');
    legend('time-aligned evt by evt', 'no time-aligned'); logy; box;
end
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (shifted & cut)'); box on;

subplot(2,2,4); hold on;
plot(tim_MWVF_cut_shifted4, PMT_MWVF_cut_shifted4); lastline('color','g');
legend('time-aligned evt by evt');
if(isPM4_StoNcut)
    plot(tim_MWVF_cut4,         PMT_MWVF_cut4);         lastline('color','g');  lastline('LineStyle','--');
    legend('time-aligned evt by evt', 'no time-aligned'); logy; box;
end
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (shifted & cut)'); box on;

figure; 
subplot(2,2,1); hold on;
plot(tim_MWVF_cut_shifted1, PMT_MWVF_cut_shifted1); lastline('color','k');
plot(tim_MWVF_cut_shifted1, 100*PMT_SWVF_cut_shifted1./PMT_MWVF_cut_shifted1); lastline('color','k');  lastline('LineStyle','--');
legend('wvf mean', 'wvf error [%]'); logy;
xlabel('time [ns]'); ylabel('amplitude [mV]'); box on;

subplot(2,2,2); hold on;
plot(tim_MWVF_cut_shifted2, PMT_MWVF_cut_shifted2); lastline('color','b');
plot(tim_MWVF_cut_shifted2, 100*PMT_SWVF_cut_shifted2./PMT_MWVF_cut_shifted2); lastline('color','b'); lastline('LineStyle','--');
legend('wvf mean', 'wvf error [%]'); logy;
xlabel('time [ns]'); ylabel('amplitude [mV]'); box on;

subplot(2,2,3); hold on;
plot(tim_MWVF_cut_shifted3, PMT_MWVF_cut_shifted3); lastline('color','r');
plot(tim_MWVF_cut_shifted3, 100*PMT_SWVF_cut_shifted3./PMT_MWVF_cut_shifted3); lastline('color','r');  lastline('LineStyle','--');
legend('wvf mean', 'wvf error [%]'); logy;
xlabel('time [ns]'); ylabel('amplitude [mV]'); box on;

subplot(2,2,4); hold on;
plot(tim_MWVF_cut_shifted4, PMT_MWVF_cut_shifted4); lastline('color','g');
plot(tim_MWVF_cut_shifted4, 100*PMT_SWVF_cut_shifted4./PMT_MWVF_cut_shifted4); lastline('color','g');  lastline('LineStyle','--');
legend('wvf mean', 'wvf error [%]'); logy;
xlabel('time [ns]'); ylabel('amplitude [mV]'); box on;


figure; hold on;
plot(tim_MWVF_cut_shifted1, PMT_MWVF_cut_shifted1/max(PMT_MWVF_cut_shifted1)); lastline('color','k');
plot(tim_MWVF_cut_shifted2, PMT_MWVF_cut_shifted2/max(PMT_MWVF_cut_shifted2)); lastline('color','b');
plot(tim_MWVF_cut_shifted3, PMT_MWVF_cut_shifted3/max(PMT_MWVF_cut_shifted3)); lastline('color','r');
plot(tim_MWVF_cut_shifted4, PMT_MWVF_cut_shifted4/max(PMT_MWVF_cut_shifted4)); lastline('color','g');
logy; yaxis (0.01,1.5);
legend('PM1', 'PM2', 'PM3', 'PM4');
xlabel('time [ns]'); ylabel('normalized amplitude '); title ('mean WVF normalized to peak (cut & shifted) '); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% VI. SINGLE-PHOTON ANALYSIS (GENERALLY NOT DONE ->BCKP!) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isSPanalysis) %#ok<BDLGI>
    A_PM1_preWVF = zeros(1,length(find(Icut))); A_PM1_postWVF = zeros(1,length(find(Icut))); A_PM1_SPcal = zeros(1,length(find(Icut)));
    A_PM2_preWVF = zeros(1,length(find(Icut))); A_PM2_postWVF = zeros(1,length(find(Icut))); A_PM2_SPcal = zeros(1,length(find(Icut)));
    A_PM3_preWVF = zeros(1,length(find(Icut))); A_PM3_postWVF = zeros(1,length(find(Icut))); A_PM3_SPcal = zeros(1,length(find(Icut)));
    A_PM4_preWVF = zeros(1,length(find(Icut))); A_PM4_postWVF = zeros(1,length(find(Icut))); A_PM4_SPcal = zeros(1,length(find(Icut)));
    
    Q_PM1_preWVF = zeros(1,length(find(Icut))); Q_PM1_postWVF = zeros(1,length(find(Icut))); Q_PM1_SPcal = zeros(1,length(find(Icut)));
    Q_PM2_preWVF = zeros(1,length(find(Icut))); Q_PM2_postWVF = zeros(1,length(find(Icut))); Q_PM2_SPcal = zeros(1,length(find(Icut)));
    Q_PM3_preWVF = zeros(1,length(find(Icut))); Q_PM3_postWVF = zeros(1,length(find(Icut))); Q_PM3_SPcal = zeros(1,length(find(Icut)));
    Q_PM4_preWVF = zeros(1,length(find(Icut))); Q_PM4_postWVF = zeros(1,length(find(Icut))); Q_PM4_SPcal = zeros(1,length(find(Icut)));
    
    for k=1:NPMTs
        icount = 0;
        for i=1:Nevts
            F = data_(k,:,i);  t = time;
            if(k==1 && Icut(i))
                icount = icount + 1;
                A_PM1_preWVF(icount)  = max(F(t<DeltaT_PM1));          Q_PM1_preWVF(icount)  = sum(F(t<DeltaT_PM1));
                A_PM1_postWVF(icount) = max(F(t>(max(t)-DeltaT_PM1))); Q_PM1_postWVF(icount) = sum(F(t>(max(t)-DeltaT_PM1)));
                A_PM1_SPcal(icount)   = max(F(t>(trigTime+Toff_PM1-DeltaT_PM1/2) & t<(trigTime+Toff_PM1+DeltaT_PM1/2))); Q_PM1_SPcal(icount)   = sum(F(t>(trigTime+Toff_PM1-DeltaT_PM1/2) & t<(trigTime+Toff_PM1+DeltaT_PM1/2)));
            end
            if(k==2 && Icut(i))
                icount = icount + 1;
                A_PM2_preWVF(icount)  = max(F(t<DeltaT_PM2));          Q_PM2_preWVF(icount)  = sum(F(t<DeltaT_PM2));
                A_PM2_postWVF(icount) = max(F(t>(max(t)-DeltaT_PM2))); Q_PM2_postWVF(icount) = sum(F(t>(max(t)-DeltaT_PM2)));
                A_PM2_SPcal(icount)   = max(F(t>(trigTime+Toff_PM2-DeltaT_PM2/2) & t<(trigTime+Toff_PM2+DeltaT_PM2/2))); Q_PM2_SPcal(icount)   = sum(F(t>(trigTime+Toff_PM2-DeltaT_PM2/2) & t<(trigTime+Toff_PM2+DeltaT_PM2/2)));
            end
            if(k==3 && Icut(i))
                icount = icount + 1;
                A_PM3_preWVF(icount)  = max(F(t<DeltaT_PM3));          Q_PM3_preWVF(icount)  = sum(F(t<DeltaT_PM3));
                A_PM3_postWVF(icount) = max(F(t>(max(t)-DeltaT_PM3))); Q_PM3_postWVF(icount) = sum(F(t>(max(t)-DeltaT_PM3)));
                A_PM3_SPcal(icount)   = max(F(t>(trigTime+Toff_PM3-DeltaT_PM3/2) & t<(trigTime+Toff_PM3+DeltaT_PM3/2))); Q_PM3_SPcal(icount)   = sum(F(t>(trigTime+Toff_PM3-DeltaT_PM3/2) & t<(trigTime+Toff_PM3+DeltaT_PM3/2)));
            end
            if(k==4 && Icut(i))
                icount = icount + 1;
                A_PM4_preWVF(icount)  = max(F(t<DeltaT_PM4));          Q_PM4_preWVF(icount)  = sum(F(t<DeltaT_PM4));
                A_PM4_postWVF(icount) = max(F(t>(max(t)-DeltaT_PM4))); Q_PM4_postWVF(icount) = sum(F(t>(max(t)-DeltaT_PM4)));
                A_PM4_SPcal(icount)   = max(F(t>(trigTime+Toff_PM4-DeltaT_PM4/2) & t<(trigTime+Toff_PM4+DeltaT_PM4/2))); Q_PM4_SPcal(icount)   = sum(F(t>(trigTime+Toff_PM4-DeltaT_PM4/2) & t<(trigTime+Toff_PM4+DeltaT_PM4/2)));
            end
        end
    end
    
    figure;
    subplot(2,2,1);
    hist1D(A_PM1_preWVF,  0:0.1:10);  lastline('color', 'k'); hold on;
    hist1D(A_PM1_postWVF, 0:0.1:10);  lastline('color', 'k'); lastline('LineStyle','--');
    xlabel('A[mV]'); legend('pre','post'); title('pre/postWVF PM1');
    subplot(2,2,2);
    hist1D(A_PM2_preWVF,  0:0.1:10);  lastline('color', 'r');  hold on;
    hist1D(A_PM2_postWVF, 0:0.1:10);  lastline('color', 'r'); lastline('LineStyle','--');
    xlabel('A[mV]'); title('pre/postWVF PM2');
    subplot(2,2,3);
    hist1D(A_PM3_preWVF,  0:0.1:10);  lastline('color', 'g');  hold on;
    hist1D(A_PM3_postWVF, 0:0.1:10);  lastline('color', 'g'); lastline('LineStyle','--');
    xlabel('A[mV]'); title('pre/postWVF PM3');
    subplot(2,2,4);
    hist1D(A_PM4_preWVF,  0:0.1:10); lastline('color', 'b');  hold on;
    hist1D(A_PM4_postWVF, 0:0.1:10);  lastline('color', 'b'); lastline('LineStyle','--');
    xlabel('A[mV]'); title('pre/postWVF PM4');
    
    figure;
    subplot(2,2,1);
    hist1D(Q_PM1_preWVF,  0:0.5:10);  lastline('color', 'k'); hold on;
    hist1D(Q_PM1_postWVF, 0:0.5:10);  lastline('color', 'k'); lastline('LineStyle','--');
    xlabel('sumA[mV]'); legend('pre','post'); title('pre/postWVF PM1');
    subplot(2,2,2);
    hist1D(Q_PM2_preWVF,  0:0.5:10);  lastline('color', 'r');  hold on;
    hist1D(Q_PM2_postWVF, 0:0.51:10);  lastline('color', 'r'); lastline('LineStyle','--');
    xlabel('sumA[mV]'); title('pre/postWVF PM2(Q)');
    subplot(2,2,3);
    hist1D(Q_PM3_preWVF,  0:0.5:10);  lastline('color', 'g');  hold on;
    hist1D(Q_PM3_postWVF, 0:0.5:10);  lastline('color', 'g'); lastline('LineStyle','--');
    xlabel('sumA[mV]'); title('pre/postWVF PM3(Q)');
    subplot(2,2,4);
    hist1D(Q_PM4_preWVF,  0:0.5:10); lastline('color', 'b');  hold on;
    hist1D(Q_PM4_postWVF, 0:0.5:10);  lastline('color', 'b'); lastline('LineStyle','--');
    xlabel('sumA[mV]'); title('pre/postWVF PM4(Q)');
    
    figure;
    subplot(2,2,1);
    hist1D(A_PM1_SPcal,  0:0.1:50);  lastline('color', 'k');
    xlabel('A[mV]'); title('SP PM1(A)');
    subplot(2,2,2);
    hist1D(A_PM2_SPcal,  0:0.1:50);  lastline('color', 'r');
    xlabel('A[mV]'); title('SP PM2(A)');
    subplot(2,2,3);
    hist1D(A_PM3_SPcal,  0:0.1:50);  lastline('color', 'g');
    xlabel('A[mV]'); title('SP PM3(A)');
    subplot(2,2,4);
    hist1D(A_PM4_SPcal,  0:0.1:50); lastline('color', 'b');
    xlabel('A[mV]'); title('SP PM4(A)');
    
    figure;
    subplot(2,2,1);
    hist1D(Q_PM1_SPcal,  0:0.5:100);  lastline('color', 'k');
    xlabel('sumA[mV]'); title('SP PM1(Q)');
    subplot(2,2,2);
    hist1D(Q_PM2_SPcal,  0:0.5:100);  lastline('color', 'r');
    xlabel('sumA[mV]'); title('SP PM2(Q)');
    subplot(2,2,3);
    hist1D(Q_PM3_SPcal,  0:0.5:100);  lastline('color', 'g');
    xlabel('sumA[mV]'); title('SP PM3(Q)');
    subplot(2,2,4);
    hist1D(Q_PM4_SPcal,  0:0.5:100); lastline('color', 'b');
    xlabel('sumA[mV]'); title('SP PM4(Q)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VII. FIT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global Ndata Tdata sigmaNdata;

figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndata      = PMT_MWVF_cut_shifted1;
sigmaNdata = PMT_SWVF_cut_shifted1;
Tdata      = tim_MWVF_cut_shifted1;
Ndata      = Ndata(Tdata<Tmax1);
sigmaNdata = sigmaNdata(Tdata<Tmax1);                                         %#ok<*NASGU>
Tdata      = Tdata(Tdata<Tmax1);                                                                                  

par0(1)   = 10;         %Amplitude of first exponential
par0(2)   = 10;         %Time constant of first exponential
par0(3)   = 10;         %Amplitude of second exponential
par0(4)   = 100;        %Time constant of second exponential
par0(5)   = 2;          %Gaussian width of PM response function
par0(6)   = 5;          %Start of exponentials (time offset)

lb(1)     = 0.1;   ub(1)     = 10000;
lb(2)     = 1;     ub(2)     = 40;
lb(3)     = 0.1;   ub(3)     = 2000;
lb(4)     = 1;     ub(4)     = 200;
lb(5)     = 1;     ub(5)     = 5;
lb(6)     = 0.3;   ub(6)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);

subplot(2,2,1); hold on;
plot(Tdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Trange = 0:0.001:1000;
y=spectrum_Exp2Conv(parT_PM1, Trange, 'PM1'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('Amplitude[mV]');
logy; yaxis(min(abs(Ndata))*0.8,max(Ndata)*1.2); xaxis(0,Tmax1);
Q_PM1_fit        = sum(spectrum_Exp2Conv(parT_PM1, Trange))*(Trange(2)-Trange(1))/Rin;                                             %charge in [pC]
Qsinglet_PM1_fit = sum(spectrum_Exp1Conv([parT_PM1(1), parT_PM1(2), parT_PM1(5), parT_PM1(6)], Trange))*(Trange(2)-Trange(1))/Rin;
Qtriplet_PM1_fit = sum(spectrum_Exp1Conv([parT_PM1(3), parT_PM1(4), parT_PM1(5), parT_PM1(6)], Trange))*(Trange(2)-Trange(1))/Rin;
Ratio_1to3_PM1   = Qsinglet_PM1_fit/Qtriplet_PM1_fit;

TmaxWVF1 = Trange(find(y==max(y)));

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndata      = PMT_MWVF_cut_shifted2;
sigmaNdata = PMT_SWVF_cut_shifted2;
Tdata      = tim_MWVF_cut_shifted2;
Ndata      = Ndata(Tdata<Tmax2);
sigmaNdata = sigmaNdata(Tdata<Tmax2);              
Tdata      = Tdata(Tdata<Tmax2);                                

par0(1)   = 10;         %Amplitude
par0(2)   = 8.5;        %Time constant of exponential
par0(3)   = 2;          %Gaussian width of PM response function
par0(4)   = 5;          %Start of exponential (time offset)

lb(1)     = 0.1;   ub(1)     = 1000;
lb(2)     = 1;     ub(2)     = 40;
lb(3)     = 1;     ub(3)     = 5;
lb(4)     = 0.3;   ub(4)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM2 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);

subplot(2,2,2); hold on;
plot(Tdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'r');
Trange = 0:0.001:1000;
y=spectrum_Exp1Conv(parT_PM2, Trange, 'PM2'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('Amplitude[mV]');
logy; yaxis(min(abs(Ndata))*0.8,max(Ndata)*1.2); xaxis(0,Tmax2);
Q_PM2_fit = sum(spectrum_Exp1Conv(parT_PM2, Trange))*(Trange(2)-Trange(1))/Rin; %charge in [pC]                

TmaxWVF2 = Trange(find(y==max(y)));

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isPM3Fit)
    Ndata      = PMT_MWVF_cut_shifted3;
    sigmaNdata = PMT_SWVF_cut_shifted3;
    Tdata      = tim_MWVF_cut_shifted3;
    Ndata      = Ndata(Tdata<Tmax3);
    sigmaNdata = sigmaNdata(Tdata<Tmax3);
    Tdata      = Tdata(Tdata<Tmax3);
    
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 5;          %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0.3;   ub(4)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    
    subplot(2,2,3); hold on;
    plot(Tdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'g');
    Trange = 0:0.001:1000;
    y=spectrum_Exp1Conv(parT_PM3, Trange, 'PM3'); xaxis(min(Trange),max(Trange));
    xlabel('time[ns]'); ylabel('Amplitude[mV]');
    logy; yaxis(min(abs(Ndata))*0.8,max(Ndata)*1.2); xaxis(0,Tmax3);
    Q_PM3_fit = sum(spectrum_Exp1Conv(parT_PM3, Trange))*(Trange(2)-Trange(1))/Rin; %charge in [pC]
    
    TmaxWVF3 = Trange(find(y==max(y)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndata      = PMT_MWVF_cut_shifted4;
sigmaNdata = PMT_SWVF_cut_shifted4;
Tdata      = tim_MWVF_cut_shifted4;
Ndata      = Ndata(Tdata<Tmax4);
sigmaNdata = sigmaNdata(Tdata<Tmax4);              
Tdata      = Tdata(Tdata<Tmax4);                                     

par0(1)   = 10;         %Amplitude of first exponential
par0(2)   = 10;         %Time constant of first exponential
par0(3)   = 10;         %Amplitude of second exponential
par0(4)   = 100;        %Time constant of second exponential
par0(5)   = 2;          %Gaussian width of PM response function
par0(6)   = 5;          %Start of exponentials (time offset)

lb(1)     = 0.1;   ub(1)     = 10000;
lb(2)     = 1;     ub(2)     = 20;
lb(3)     = 0.1;   ub(3)     = 2000;
lb(4)     = 1;     ub(4)     = 200;
lb(5)     = 1;     ub(5)     = 5;
lb(6)     = 0.3;   ub(6)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM4 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);

subplot(2,2,4); hold on;
plot(Tdata, Ndata, '+', 'Markersize', 10);  lastline('color', 'k');
Trange = 0:0.001:1000;
y=spectrum_Exp2Conv(parT_PM4, Trange, 'PM4'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('Amplitude[mV]');
logy; yaxis(min(abs(Ndata))*0.8,max(Ndata)*1.2); xaxis(0,Tmax4);
Ratio_1to3_PM4   = sum(spectrum_Exp1Conv([parT_PM4(1), parT_PM4(2), parT_PM4(5), parT_PM4(6)], Trange))/sum(spectrum_Exp1Conv([parT_PM4(3), parT_PM4(4), parT_PM4(5), parT_PM4(6)], Trange));
Q_PM4_fit        = sum(spectrum_Exp2Conv(parT_PM4, Trange))*(Trange(2)-Trange(1))/Rin; %charge in [pC]
Qsinglet_PM4_fit = sum(spectrum_Exp1Conv([parT_PM4(1), parT_PM4(2), parT_PM4(5), parT_PM4(6)], Trange))*(Trange(2)-Trange(1))/Rin; %charge in [pC]
Qtriplet_PM4_fit = sum(spectrum_Exp1Conv([parT_PM4(3), parT_PM4(4), parT_PM4(5), parT_PM4(6)], Trange))*(Trange(2)-Trange(1))/Rin; %charge in [pC]
Ratio_1to3_PM4   = Qsinglet_PM4_fit/Qtriplet_PM4_fit;

TmaxWVF4 = Trange(find(y==max(y)));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% VIIb. OBTAIN CALIBRATED DISTRIBUTIONS AND ALIGN: %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data___(1,:,:) = (data__(1,:,:)*TBin/Rin)/QPM1toPh*1/TBin;                  %phe/ns
data___(2,:,:) = (data__(2,:,:)*TBin/Rin)/QPM2toPh*1/TBin;                  %phe/ns
data___(3,:,:) = (data__(3,:,:)*TBin/Rin)/QPM3toPh*1/TBin;                  %phe/ns
data___(4,:,:) = (data__(4,:,:)*TBin/Rin)/QPM4toPh*1/TBin;                  %phe/ns

wvf1  = (PMT_MWVF_cut_shifted1*TBin/Rin)/QPM1toPh*1/TBin;                   %phe/ns
wvf2  = (PMT_MWVF_cut_shifted2*TBin/Rin)/QPM2toPh*1/TBin;                   %phe/ns
wvf3  = (PMT_MWVF_cut_shifted3*TBin/Rin)/QPM3toPh*1/TBin;                   %phe/ns
wvf4  = (PMT_MWVF_cut_shifted4*TBin/Rin)/QPM4toPh*1/TBin;                   %phe/ns

swvf1 = (PMT_SWVF_cut_shifted1*TBin/Rin)/QPM1toPh*1/TBin;                   %phe/ns
swvf2 = (PMT_SWVF_cut_shifted2*TBin/Rin)/QPM2toPh*1/TBin;                   %phe/ns
swvf3 = (PMT_SWVF_cut_shifted3*TBin/Rin)/QPM3toPh*1/TBin;                   %phe/ns
swvf4 = (PMT_SWVF_cut_shifted4*TBin/Rin)/QPM4toPh*1/TBin;                   %phe/ns

wvf1  = interp1(tim_MWVF_cut_shifted1,  wvf1, time, 'linear', 0);
swvf1 = interp1(tim_MWVF_cut_shifted1, swvf1, time, 'linear', 0);
tf1   = TmaxWVF1; ti1 = TmaxWVF1;
wvf1  = TimeAlign(time,  wvf1, ti1, tf1);
swvf1 = TimeAlign(time, swvf1, ti1, tf1);

wvf2  = interp1(tim_MWVF_cut_shifted2,  wvf2, time, 'linear', 0);
swvf2 = interp1(tim_MWVF_cut_shifted2, swvf2, time, 'linear', 0);
tf2    = TmaxWVF1; ti2 = TmaxWVF2;
wvf2  = TimeAlign(time,  wvf2, ti2, tf2);
swvf2 = TimeAlign(time, swvf2, ti2, tf2);

wvf3  = interp1(tim_MWVF_cut_shifted3,  wvf3, time, 'linear', 0);
swvf3 = interp1(tim_MWVF_cut_shifted3, swvf3, time, 'linear', 0);
if(isPM3Fit),    tf3    = TmaxWVF1; ti3 = TmaxWVF3;
else,            tf3    = TmaxWVF1; ti3 = time(min(find(wvf3==max(wvf3))));
end
wvf3  = TimeAlign(time,  wvf3, ti3, tf3);
swvf3 = TimeAlign(time, swvf3, ti3, tf3);

wvf4  = interp1(tim_MWVF_cut_shifted4,  wvf4, time, 'linear', 0);
swvf4 = interp1(tim_MWVF_cut_shifted4, swvf4, time, 'linear', 0);
tf4    = TmaxWVF1; ti4 = TmaxWVF4;
wvf4  = TimeAlign(time,  wvf4, ti4, tf4);
swvf4 = TimeAlign(time, swvf4, ti4, tf4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% VIIc. FIT, REMOVE AND (OPTIONALLY) STORE AFTERPULSING TEMPLATE %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvf1_  = (PMT_MWVF_cut_shifted(1,:)*TBin/Rin)/QPM1toPh*1/TBin;              %phe/ns
wvf2_  = (PMT_MWVF_cut_shifted(2,:)*TBin/Rin)/QPM2toPh*1/TBin;              %phe/ns
wvf3_  = (PMT_MWVF_cut_shifted(3,:)*TBin/Rin)/QPM3toPh*1/TBin;              %phe/ns
wvf4_  = (PMT_MWVF_cut_shifted(4,:)*TBin/Rin)/QPM4toPh*1/TBin;              %phe/ns

swvf1_ = (PMT_SWVF_cut_shifted(1,:)*TBin/Rin)/QPM1toPh*1/TBin;              %phe/ns
swvf2_ = (PMT_SWVF_cut_shifted(2,:)*TBin/Rin)/QPM2toPh*1/TBin;              %phe/ns
swvf3_ = (PMT_SWVF_cut_shifted(3,:)*TBin/Rin)/QPM3toPh*1/TBin;              %phe/ns
swvf4_ = (PMT_SWVF_cut_shifted(4,:)*TBin/Rin)/QPM4toPh*1/TBin;              %phe/ns

if(storeAfterPulsingTemplate)   %Do now fit to calibrated data to obtain afterpulsing template at times above 200ns
    
    figure;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PMT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,1);
    Ndata      = wvf1_;
    sigmaNdata = swvf1_;
    sigmaNdata(find(sigmaNdata==0)) = max(sigmaNdata);                      %Reassign errors to zeros
    Tdata      = time;
    Ndata      = Ndata(Tdata<Tmax1);
    sigmaNdata = sigmaNdata(Tdata<Tmax1);                                   %#ok<*NASGU>
    Tdata      = Tdata(Tdata<Tmax1);
    
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 5;          %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 0.3;   ub(6)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_tmp = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    y=spectrum_Exp2Conv(parT_tmp, time); xaxis(min(Trange),max(Trange));
    afterp1 = (wvf1_-y)/max(wvf1_);
    afterp1(time<=190) = 0;
    plot(time, afterp1, '-', 'Markersize', 10);  lastline('color', 'k'); box; box;
    xlabel('time[ns]'); ylabel('afterpulsing (normalized to maximum)');
    yaxis(max(afterp1)/100,max(afterp1)*1.2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PMT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,2);
    Ndata      = wvf2_;
    sigmaNdata = swvf2_;
    sigmaNdata(find(sigmaNdata==0)) = max(sigmaNdata);                      %Reassign errors to zeros
    Tdata      = time;
    Ndata      = Ndata(Tdata<Tmax2);
    sigmaNdata = sigmaNdata(Tdata<Tmax2);                                   %#ok<*NASGU>
    Tdata      = Tdata(Tdata<Tmax2);
    
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 5;          %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0.3;   ub(4)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_tmp = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    y=spectrum_Exp1Conv(parT_tmp, time);
    afterp2 = (wvf2_-y)/max(wvf2_);
    afterp2(time<=190) = 0;
    plot(time, afterp2, '-', 'Markersize', 10);  lastline('color', 'k'); box; box;
    xlabel('time[ns]'); ylabel('afterpulsing (normalized to maximum)');
    yaxis(max(afterp2)/100,max(afterp2)*1.2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PMT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,3);
    Ndata      = wvf3_;
    sigmaNdata = swvf3_;
    sigmaNdata(find(sigmaNdata==0)) = max(sigmaNdata);                      %Reassign errors to zeros
    Tdata      = time;
    Ndata      = Ndata(Tdata<Tmax3);
    sigmaNdata = sigmaNdata(Tdata<Tmax3);                                   %#ok<*NASGU>
    Tdata      = Tdata(Tdata<Tmax3);
    
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 5;          %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0.3;   ub(4)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_tmp = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    y=spectrum_Exp1Conv(parT_tmp, time);
    afterp3 = (wvf3_-y)/max(wvf3_);
    afterp3(time<=190) = 0;
    plot(time, afterp3, '-', 'Markersize', 10);  lastline('color', 'k'); box; box;
    xlabel('time[ns]'); ylabel('afterpulsing (normalized to maximum)');
    yaxis(max(afterp3)/100,max(afterp3)*1.2);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%% PMT4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,2,4);
    Ndata      = wvf4_;
    sigmaNdata = swvf4_;
    sigmaNdata(find(sigmaNdata==0)) = max(sigmaNdata);                      %Reassign errors to zeros
    Tdata      = time;
    Ndata      = Ndata(Tdata<Tmax4);
    sigmaNdata = sigmaNdata(Tdata<Tmax4);                                   %#ok<*NASGU>
    Tdata      = Tdata(Tdata<Tmax4);
    
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 5;          %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 0.3;   ub(6)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_tmp = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    y=spectrum_Exp2Conv(parT_tmp, time); xaxis(min(Trange),max(Trange));
    afterp4 = (wvf4_-y)/max(wvf4_);
    afterp4(time<=190) = 0;
    plot(time, afterp4, '-', 'Markersize', 10);  lastline('color', 'k'); box; box;
    xlabel('time[ns]'); ylabel('afterpulsing (normalized to maximum)');
    yaxis(max(afterp4)/100,max(afterp4)*1.2);
end

if(storeAfterPulsingTemplate)   %Climb-up to parent folder. FIXME: it would be more general if storing the time_shift.    
    i=1;
    while ((DIR_OUTPUT(length(DIR_OUTPUT)-i)) ~= '\'), i=i+1; end
    DIR_OUTPUTa = DIR_OUTPUT(1:(length(DIR_OUTPUT)-i));
    fid  = fopen([DIR_OUTPUTa, 'afterpulse'], 'w');
    fprintf(fid, '%s\r', 'time afterp1 afterp2 afterp3 afterp4');
    for j=1:length(time)
        fprintf(fid, '%3.7f %3.7f %3.7f %3.7f %3.7f \r', time(j), afterp1(j), afterp2(j), afterp3(j), afterp4(j));
    end
    fclose('all');
end

if(useAfterPulsingTemplate)
    if(~TakeAfterPulsingFromExt), DIR_OUTPUTa = DIR_OUTPUT(1:70); end       %FIXME: not sure how this is suppossed to work
    
    fid  = fopen([DIR_OUTPUTa, 'afterpulse'], 'r'); 
    fgets(fid);
    A = fscanf(fid, '%f %f %f %f %f',[5, inf]);
    time    = A(1,:);
    afterp1 = A(2,:);
    afterp2 = A(3,:);    
    afterp3 = A(4,:);    
    afterp4 = A(5,:);
    fclose('all');
else
    afterp1 = zeros(size(time));
    afterp2 = zeros(size(time)); 
    afterp3 = zeros(size(time));   
    afterp4 = zeros(size(time));
end

if(useAfterPulsingTemplate)
figure;
subplot(2,2,1); 
plot(time, wvf1_,'-'); logy; yaxis(0.001*max(wvf1_), 1.2*max(wvf1_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf with afterpulsing');
subplot(2,2,2); 
plot(time, wvf2_,'-'); logy; yaxis(0.001*max(wvf2_), 1.2*max(wvf2_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf with afterpulsing');
subplot(2,2,3); 
plot(time, wvf3_,'-'); logy; yaxis(0.001*max(wvf3_), 1.2*max(wvf3_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf with afterpulsing');
subplot(2,2,4); 
plot(time, wvf4_,'-'); logy; yaxis(0.001*max(wvf4_), 1.2*max(wvf4_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf with afterpulsing');
    
figure;
subplot(2,2,1); 
plot(time, (wvf1_ - afterp1*max(wvf1_)),'-'); logy; yaxis(0.001*max(wvf1_), 1.2*max(wvf1_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf without afterpulsing');
subplot(2,2,2); 
plot(time, (wvf2_ - afterp2*max(wvf2_)),'-'); logy; yaxis(0.001*max(wvf2_), 1.2*max(wvf2_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf without afterpulsing');
subplot(2,2,3); 
plot(time, (wvf3_ - afterp3*max(wvf3_)),'-'); logy; yaxis(0.001*max(wvf3_), 1.2*max(wvf3_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf without afterpulsing');
subplot(2,2,4); 
plot(time, (wvf4_ - afterp4*max(wvf4_)),'-'); logy; yaxis(0.001*max(wvf4_), 1.2*max(wvf4_));
ylabel('phe/ns'); xlabel('time[ns]'); title('wvf without afterpulsing');

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VIId. FIT FINAL WVFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvf1__ = wvf1_ - afterp1*max(wvf1_);
wvf2__ = wvf2_ - afterp2*max(wvf2_);
wvf3__ = wvf3_ - afterp3*max(wvf3_);
wvf4__ = wvf4_ - afterp4*max(wvf4_);

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

Ndata      = wvf1__;
sigmaNdata = swvf1_;
Tdata      = time;
Ndata      = Ndata(Tdata<Tmax1);
sigmaNdata = sigmaNdata(Tdata<Tmax1);                                         %#ok<*NASGU>
Tdata      = Tdata(Tdata<Tmax1);                                                                                  

par0(1)   = 10;         %Amplitude of first exponential
par0(2)   = 10;         %Time constant of first exponential
par0(3)   = 10;         %Amplitude of second exponential
par0(4)   = 100;        %Time constant of second exponential
par0(5)   = 2;          %Gaussian width of PM response function
par0(6)   = 5;          %Start of exponentials (time offset)

lb(1)     = 0.1;   ub(1)     = 10000;
lb(2)     = 1;     ub(2)     = 40;
lb(3)     = 0.1;   ub(3)     = 2000;
lb(4)     = 1;     ub(4)     = 200;
lb(5)     = 1;     ub(5)     = 5;
lb(6)     = 0.3;   ub(6)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM1_ = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);

subplot(2,2,1); hold on;
plot(time, wvf1__, '+', 'Markersize', 10);  lastline('color', 'k');
Trange = 0:0.001:1000;
y=spectrum_Exp2Conv(parT_PM1_, Trange, 'PM1 (final)'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('spectrum [pe/ns]');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(1,4500);logy; logx; 

Q_PM1_fit_pe_        = sum(spectrum_Exp2Conv(parT_PM1_, Trange))*(Trange(2)-Trange(1));
Qsinglet_PM1_fit_pe_ = sum(spectrum_Exp1Conv([parT_PM1_(1), parT_PM1_(2), parT_PM1_(5), parT_PM1_(6)], Trange))*(Trange(2)-Trange(1));
Qtriplet_PM1_fit_pe_ = sum(spectrum_Exp1Conv([parT_PM1_(3), parT_PM1_(4), parT_PM1_(5), parT_PM1_(6)], Trange))*(Trange(2)-Trange(1));
Ratio_1to3_PM1_pe_   = Qsinglet_PM1_fit_pe_/Qtriplet_PM1_fit_pe_;

Q_PM1_all_pe_ = sum(wvf1__)*TBin;
max1tmp = max(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndata      = wvf2__;
sigmaNdata = swvf2_;
Tdata      = time;
Ndata      = Ndata(Tdata<Tmax2);
sigmaNdata = sigmaNdata(Tdata<Tmax2);              
Tdata      = Tdata(Tdata<Tmax2);                                

par0(1)   = 10;         %Amplitude
par0(2)   = 8.5;        %Time constant of exponential
par0(3)   = 2;          %Gaussian width of PM response function
par0(4)   = 5;          %Start of exponential (time offset)

lb(1)     = 0.1;   ub(1)     = 1000;
lb(2)     = 1;     ub(2)     = 40;
lb(3)     = 1;     ub(3)     = 5;
lb(4)     = 0.3;   ub(4)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM2_ = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);

subplot(2,2,2); hold on;
plot(time, wvf2__, '+', 'Markersize', 10); lastline('color', 'r');
Trange = 0:0.001:1000;
y=spectrum_Exp1Conv(parT_PM2_, Trange, 'PM2 (final)'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('spectrum [pe/ns]');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(1,4500);logy; logx; 

Q_PM2_fit_pe_ = sum(spectrum_Exp1Conv(parT_PM2_, Trange))*(Trange(2)-Trange(1));               

Q_PM2_Tcut_pe_  = sum(wvf2__(time<Tmax2))*TBin;
max2tmp = max(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(isPM3Fit)
    Ndata      = wvf3__;
    sigmaNdata = swvf3_;
    Tdata      = time;
    Ndata      = Ndata(Tdata<Tmax3);
    sigmaNdata = sigmaNdata(Tdata<Tmax3);
    Tdata      = Tdata(Tdata<Tmax3);
    
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 5;          %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0.3;   ub(4)     = 15;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3_ = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    
    subplot(2,2,3); hold on;
    plot(time, wvf3__, '+', 'Markersize', 10);  lastline('color', 'g');
    Trange = 0:0.001:1000;
    y=spectrum_Exp1Conv(parT_PM3_, Trange, 'PM3 (final)'); xaxis(min(Trange),max(Trange));
    xlabel('time[ns]'); ylabel('spectrum [pe/ns]');
    yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(1,4500);logy; logx;
    
    Q_PM3_fit_pe_ = sum(spectrum_Exp1Conv(parT_PM3_, Trange))*(Trange(2)-Trange(1));
    max3tmp = max(y);
end
Q_PM3_Tcut_pe_  = sum(wvf3__(time<Tmax3))*TBin;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ndata      = wvf4__;
sigmaNdata = swvf4_;
Tdata      = time;
Ndata      = Ndata(Tdata<Tmax4);
sigmaNdata = sigmaNdata(Tdata<Tmax4);              
Tdata      = Tdata(Tdata<Tmax4);                                     

par0(1)   = 10;         %Amplitude of first exponential
par0(2)   = 10;         %Time constant of first exponential
par0(3)   = 10;         %Amplitude of second exponential
par0(4)   = 1000;       %Time constant of second exponential
par0(5)   = 2;          %Gaussian width of PM response function
par0(6)   = 5;          %Start of exponentials (time offset)

%lb(1)     = 0.1;   ub(1)     = 10000;
%lb(2)     = 0.1;   ub(2)     = 20;
%lb(3)     = 0.1;   ub(3)     = 5000; %PERHAPS ADJUST FOR PURE AR
%lb(4)     = 400;   ub(4)     = 3000;
%lb(5)     = 1;     ub(5)     = 5;
%lb(6)     = 0.3;   ub(6)     = 15;

lb(1)     = 0.1;   ub(1)     = 10000;
lb(2)     = 1;     ub(2)     = 40;
lb(3)     = 0.1;   ub(3)     = 2000;
lb(4)     = 1;     ub(4)     = 200;
lb(5)     = 1;     ub(5)     = 5;
lb(6)     = 0.3;   ub(6)     = 15;

options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
parT_PM4_ = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);

subplot(2,2,4); hold on;
plot(time, wvf4__, '+', 'Markersize', 10);  lastline('color', 'k');
Trange = 0:0.001:1000;
y=spectrum_Exp2Conv(parT_PM4_, Trange, 'PM4 (final)'); xaxis(min(Trange),max(Trange));
xlabel('time[ns]'); ylabel('spectrum [pe/ns]');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(1,4500);logy; logx; 

Ratio_1to3_PM4_      = sum(spectrum_Exp1Conv([parT_PM4_(1), parT_PM4_(2), parT_PM4_(5), parT_PM4_(6)], Trange))/sum(spectrum_Exp1Conv([parT_PM4_(3), parT_PM4_(4), parT_PM4_(5), parT_PM4_(6)], Trange));
Q_PM4_fit_pe_        = sum(spectrum_Exp2Conv(parT_PM4_, Trange))*(Trange(2)-Trange(1));
Qsinglet_PM4_fit_pe_ = sum(spectrum_Exp1Conv([parT_PM4_(1), parT_PM4_(2), parT_PM4_(5), parT_PM4_(6)], Trange))*(Trange(2)-Trange(1));
Qtriplet_PM4_fit_pe_ = sum(spectrum_Exp1Conv([parT_PM4_(3), parT_PM4_(4), parT_PM4_(5), parT_PM4_(6)], Trange))*(Trange(2)-Trange(1));
Ratio_1to3_PM4_pe_   = Qsinglet_PM4_fit_pe_/Qtriplet_PM4_fit_pe_;

Q_PM4_all_pe_ = sum(wvf4__)*TBin;

max4tmp = max(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% VIII. CALCULATE W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%(correct for solid angle, quantum efficiency, transparency.
% For Xe 2nd continuum: correct for quenching)

Q_PM1_Tcut          = sum(PMT_MWVF_cut_shifted1(tim_MWVF_cut_shifted1<Tmax1))*TBin/Rin;
Q_PM1_all	        = sum(PMT_MWVF_cut(1,:))*TBin/Rin;
Q_PM2_Tcut          = sum(PMT_MWVF_cut_shifted2(tim_MWVF_cut_shifted2<Tmax2))*TBin/Rin;
Q_PM2_all	        = sum(PMT_MWVF_cut(2,:))*TBin/Rin;
Q_PM3_Tcut          = sum(PMT_MWVF_cut_shifted3(tim_MWVF_cut_shifted3<Tmax3))*TBin/Rin;
Q_PM3_all	        = sum(PMT_MWVF_cut(3,:))*TBin/Rin;
Q_PM4_Tcut          = sum(PMT_MWVF_cut_shifted4(tim_MWVF_cut_shifted4<Tmax4))*TBin/Rin;
Q_PM4_all	        = sum(PMT_MWVF_cut(4,:))*TBin/Rin;

Q_PM1_fit_pe        = Q_PM1_fit/QPM1toPh;
Q_PM1_Tcut_pe       = Q_PM1_Tcut/QPM1toPh;
Q_PM1_all_pe        = Q_PM1_all/QPM1toPh;

Q_PM2_fit_pe        = Q_PM2_fit/QPM2toPh;
Q_PM2_Tcut_pe       = Q_PM2_Tcut/QPM2toPh;
Q_PM2_all_pe        = Q_PM2_all/QPM2toPh;

if(isPM3Fit), Q_PM3_fit_pe        = Q_PM3_fit/QPM3toPh; end
Q_PM3_Tcut_pe       = Q_PM3_Tcut/QPM3toPh;
Q_PM3_all_pe        = Q_PM3_all/QPM3toPh;

Q_PM4_fit_pe        = Q_PM4_fit/QPM4toPh;
Q_PM4_Tcut_pe       = Q_PM4_Tcut/QPM4toPh;
Q_PM4_all_pe        = Q_PM4_all/QPM4toPh;

Ealpha   = 5.5*1e+6;                  %Energy of alpha particle.

Omega    = 0.0085 * 1.72; %Makes average between filter and no filter (0.00875 and 0.0081)for nominal (Xe-paper) geometry. 
                          %Introduces a factor 1.72 to reproduce Xe-Wsc (simulation by Angela gives ~2.7, 
                          %but not in exact same geometry. It should be less, perhaps 2.)

if (isXe), QE_PM1   = 0.137; end       %PM7378 efficiency for Xe 2nd continuum
if (isAr), QE_PM1   = 0.185; end       %PM7378 efficiency for Ar 3rd continuum
if (~isXe && ~isAr), QE_PM1  = 0.21;   end       %PM7378 efficiency (typical -> take the same as PM2)

if (isXe), QE_PM2   = 0.182; end       %PM7378 efficiency for Xe 3rd continuum with 250-400 filter
if (isAr), QE_PM2   = 0.049; end       %PM7378 efficiency for Ar 3rd continuum with 250-400 filter
if (~isXe && ~isAr), QE_PM2   = 0.21;  end       %PM7378 efficiency (from Pablo's spectrum: 1%CF4->0.21, 10%CF4->0.2095)

if (isXe), QE_PM3   = 0.0788; end      %PM5070 efficiency for Xe 3rd continuum with 400-700 filter
if (isAr), QE_PM3   = 6.9770e-05; end  %PM5070 efficiency for Ar 3rd continuum with 400-700 filter
%OLD ASSUMPTION
%QE_PM3   = 0.1;                      %PM5070 efficiency in range 400-700,
%NEW CALCULATION
if (~isXe && ~isAr), QE_PM3   = 0.0802; end    %PM5070 efficiency (from Pablo's spectrum: 1%CF4->0.0816, 10%CF4->0.0788)

Omega_wfilter    = 0.0085 * 1.72; %Makes average between filter and no filter (0.00875 and 0.0081)for nominal (Xe-paper) geometry. 
                                  %Introduces a factor 1.72 to reproduce Xe-Wsc (simulation by Angela gives ~2.7, 
                                  %but not in exact same geometry. It should be less, perhaps 2.)
                                  
QE_400nm_PM4 = 0.25;                    %Quantum efficiency of PM4 at 400nm
TPB_WLE      = 0.2;                     %WLE of TPB
OmegaTPB     = 0.3;                     %Approximate solid angle of TPB at PM entrance
QE_PM4       = TPB_WLE * OmegaTPB * QE_400nm_PM4; %Quantum efficiency of PM4 to UV light

Transp  = 0.936*(0.968^2)*0.99;      %Transparency of anode (thin) and 3 additional planes.

nph_prod = Q_PM2_Tcut_pe_ /(Omega_wfilter*QE_PM2*Transp);
W2_Tcut = Ealpha/nph_prod;
nph_prod = Q_PM2_fit_pe_ /(Omega_wfilter*QE_PM2*Transp);
W2_fit = Ealpha/nph_prod;

nph_prod = Q_PM3_Tcut_pe_ /(Omega_wfilter*QE_PM3*Transp);
W3_Tcut = Ealpha/nph_prod;

if(isPM3Fit)
    nph_prod = Q_PM3_fit_pe_ /(Omega_wfilter*QE_PM3*Transp);
    W3_fit = Ealpha/nph_prod;
end

nph_prod = (Q_PM1_all_pe_ + isCorr*Qtriplet_PM1_fit_pe_*(100/parT_PM1_(4)-1))/(Omega*QE_PM1*Transp);
W1       = Ealpha/nph_prod;
nph_prod = (Q_PM1_fit_pe_ + isCorr*Qtriplet_PM1_fit_pe_*(100/parT_PM1_(4)-1))/(Omega*QE_PM1*Transp);
W1_fit   = Ealpha/nph_prod;

nph_prod = (Q_PM4_all_pe_ + isCorr*Qtriplet_PM4_fit_pe_*(100/parT_PM4_(4)-1))/(Omega*QE_PM4*Transp);
W4       = Ealpha/nph_prod;
nph_prod = (Q_PM4_fit_pe_ + isCorr*Qtriplet_PM4_fit_pe_*(100/parT_PM4_(4)-1))/(Omega*QE_PM4*Transp);
W4_fit   = Ealpha/nph_prod;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% IX. SUMMARY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('                                                                   ');
disp('------------------------   summary   ------------------------------');
disp('                                                                   ');
if(useAfterPulsingTemplate)
    disp('        ***Note: afterpulsing has been subtracted!***       ')
end

disp('                PM1 R7378 (no filter)                    ');

if(~useAfterPulsingTemplate)    
    disp(['nph detected (fit)          = ', num2str(Q_PM1_fit_pe),  ' (PM1)']);
    disp(['nph detected [<180ns]       = ', num2str(Q_PM1_Tcut_pe),  ' (PM1)']);
    disp(['tau1                        = ', num2str(parT_PM1(2)),  ' ns        ', 'sigma0 = ', num2str(parT_PM1(5)),  ' ns']);
    disp(['tau3                        = ', num2str(parT_PM1(4)),  ' ns        ', 't0     = ', num2str(parT_PM1(6)), ' ns']);
else
    disp(['nph detected (fit)          = ', num2str(Q_PM1_fit_pe),  ' (PM1)']);
    disp(['nph detected [<180ns]       = ', num2str(Q_PM1_Tcut_pe),  ' (PM1)']);
    disp(['tau1                        = ', num2str(parT_PM1_(2)),  ' ns       ', 'sigma0 = ', num2str(parT_PM1_(5)),  ' ns']);
    disp(['tau3                        = ', num2str(parT_PM1_(4)),  ' ns       ', 't0     = ', num2str(parT_PM1_(6)), ' ns']);
end

disp(['Wsc-value (PM1)                      = ', num2str(W1),     ' [eV] (PM1)']);
disp(['Wsc-value (PM1), fit                 = ', num2str(W1_fit), ' [eV] (PM1)']);

disp('                                                         ');
disp('         PM2 R7378 (filter: >250 nm, <400 nm)            ');

if(~useAfterPulsingTemplate)    
    disp(['nph detected (fit)          = ', num2str(Q_PM2_fit_pe),  ' (PM2)']);
    disp(['nph detected [<50ns]        = ', num2str(Q_PM2_Tcut_pe),   ' (PM2)']);
    disp(['tau                         = ', num2str(parT_PM2(2)),   ' ns       ', 'sigma0 = ', num2str(parT_PM2(5)),  ' ns']);
else
    disp(['nph detected (fit)          = ', num2str(Q_PM2_fit_pe),  ' (PM2)']);
    disp(['nph detected [<50ns]        = ', num2str(Q_PM2_Tcut_pe),   ' (PM2)']);
    disp(['tau                         = ', num2str(parT_PM2_(2)),  ' ns       ', 'sigma0 = ', num2str(parT_PM2_(5)),  ' ns']);
end

disp(['Wsc-value (PM2)                      = ', num2str(W2_Tcut),     ' [eV] (PM2)']);
disp(['Wsc-value (PM2), fit                 = ', num2str(W2_fit), ' [eV] (PM2)']);

disp('                                                         ');
disp('         PM3 R5070 (filter: >400 nm, <700 nm)            ');

if(~useAfterPulsingTemplate)    
    if(isPM3Fit), disp(['nph detected (fit)          = ', num2str(Q_PM3_fit_pe),  ' (PM1)']); end
    disp(['nph detected [<50ns]        = ', num2str(Q_PM3_Tcut_pe),   ' (PM3)']);
    if(isPM3Fit), disp(['tau                         = ', num2str(parT_PM3(2)),   ' ns       ', 'sigma0 = ', num2str(parT_PM3(5)),  ' ns']); end
else
    if(isPM3Fit), disp(['nph detected (fit)          = ', num2str(Q_PM3_fit_pe),  ' (PM1)']); end
    disp(['nph detected [<50ns]        = ', num2str(Q_PM3_Tcut_pe),   ' (PM3)']);
    if(isPM3Fit), disp(['tau                         = ', num2str(parT_PM3_(2)),  ' ns       ', 'sigma0 = ', num2str(parT_PM3_(5)),  ' ns']); end
end

disp(['Wsc-value (PM3)                      = ', num2str(W3_Tcut),     ' [eV] (PM3)']);
if(isPM3Fit), disp(['Wsc-value (PM3), fit                 = ', num2str(W3_fit), ' [eV] (PM3)']); end

disp('                                                                   ');
disp('                PM4 R7378 (with TPB)                    ');

if(~useAfterPulsingTemplate)    
    disp(['nph detected (fit)          = ', num2str(Q_PM4_fit_pe),  ' (PM4)']);
    disp(['nph detected [<200ns]       = ', num2str(Q_PM4_Tcut_pe),  ' (PM4)']);
    disp(['tau1                        = ', num2str(parT_PM4(2)),  ' ns        ', 'sigma0 = ', num2str(parT_PM4(5)),  ' ns']);
    disp(['tau3                        = ', num2str(parT_PM4(4)),  ' ns        ', 't0     = ', num2str(parT_PM4(6)), ' ns']);
else
    disp(['nph detected (fit)          = ', num2str(Q_PM4_fit_pe),  ' (PM4)']);
    disp(['nph detected [<200ns]       = ', num2str(Q_PM4_Tcut_pe),  ' (PM4)']);
    disp(['tau1                        = ', num2str(parT_PM4_(2)),  ' ns       ', 'sigma0 = ', num2str(parT_PM4_(5)),  ' ns']);
    disp(['tau3                        = ', num2str(parT_PM4_(4)),  ' ns       ', 't0     = ', num2str(parT_PM4_(6)), ' ns']);
end

disp(['Wsc-value (PM4)                      = ', num2str(W4),     ' [eV] (PM4)']);
disp(['Wsc-value (PM4), fit                 = ', num2str(W4_fit), ' [eV] (PM4)']);

disp('                                                   ');
disp(['total number of events = ', num2str(Nevts)]);
disp(['selected events        = ', num2str(length(find(Icut)))]);
disp('                                                   ');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% X. SAVE TO FILE: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(~isPM3Fit && (isSaveSummaryToFile || isSaveWVFsToFile)) %#ok<BDLGI>
    disp('****WARNING****: FILE NOT SAVED SINCE PM3 NOT FITTED!');
elseif (isPM3Fit) %#ok<BDLGI>
    if(isSaveSummaryToFile)
        %Trim waveforms (calibrated, shifted and afterpulse-suppressed) to interesting time-range:
        Tmax = max([Tmax1; Tmax2; Tmax3; Tmax4]);
        wvf1_SAVE = wvf1__(time<Tmax); swvf1_SAVE = swvf1_(time<Tmax);
        wvf2_SAVE = wvf2__(time<Tmax); swvf2_SAVE = swvf2_(time<Tmax);
        wvf3_SAVE = wvf3__(time<Tmax); swvf3_SAVE = swvf3_(time<Tmax);
        wvf4_SAVE = wvf4__(time<Tmax); swvf4_SAVE = swvf4_(time<Tmax);
        time_SAVE = time(time<Tmax);
        
        fid   = fopen([DIR_OUTPUT, FILE, '_cal'], 'w');
        fprintf(fid, '%s\r', 'time   N1[ped/ns]   N2[ped/ns]   N3[ped/ns]   N4[ped/ns]   s_N1[ped/ns]   s_ N2[ped/ns]   s_N3[ped/ns]   s_N4[ped/ns]');
        
        for j=1:length(time_SAVE)
            fprintf(fid, '%3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f \r', time_SAVE(j), wvf1_SAVE(j), wvf2_SAVE(j), wvf3_SAVE(j), ...
                wvf4_SAVE(j), swvf1_SAVE(j), swvf2_SAVE(j), swvf3_SAVE(j), swvf4_SAVE(j));
        end
        
        %            %THIS MIGHT BE GOOD TO SAVE SUMMARY PARAMETERS
        %            fid2   = fopen([DIR_OUTPUT, FILE, '_data'], 'w');
        %            fprintf(fid2, '%3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f, %3.3f',...
        %            Q_PM1_fit_pe, Q_PM1_fit_pe_corr, Q_PM1_all_pe_, parT_PM1(2), parT_PM1(4),...
        %            Q_PM4_fit_pe, Q_PM4_fit_pe_corr, Q_PM4_all_pe_, parT_PM4(2), parT_PM4(4),...
        %            Q_PM2_fit_pe, Q_PM2_Tcut_pe, parT_PM2(2), Q_PM3_fit_pe, Q_PM3_Tcut_pe, parT_PM3(2),...
        %            W1/W2_Tcut, W1/W3_Tcut,...
        %            W1, W4, Nevts, length(find(Icut)));
    end
    if (isSaveWVFsToFile) %#ok<BDLGI>
        
        %trim waveforms (calibrated, shifted and afterpulse-suppressed) to interesting time-range:
        Tmax = max([Tmax1; Tmax2; Tmax3; Tmax4]);
        
        fid3  = fopen([DIR_OUTPUT, FILE, '_wvfs_cal'], 'w');
        fprintf(fid3, '%s\r', 'time   N1[ped/ns]   N2[ped/ns]   N3[ped/ns]   N4[ped/ns] ');
        
        Nevts_cut = length(A1(Icut));
        if(NwvfsWritten && NwvfsWritten < Nevts_cut), Nevts_cut = NwvfsWritten; end
        for i=1:Nevts_cut
            fprintf(fid3, '%s\r', num2str(i));
            if(rem(i, 100)==0), disp(['wvfs ', num2str(i), ' written to disk']); end
            for j=1:length(time)
                fprintf(fid3, '%3.3f, %3.3f, %3.3f, %3.3f, %3.3f \r', time(j), data___(1,j,i), data___(2,j,i), data___(3,j,i), data___(4,j,i));
                if (time(j)>Tmax), break; end
            end
        end
    end
end

fclose('all');
toc;

mosaic;
return;

