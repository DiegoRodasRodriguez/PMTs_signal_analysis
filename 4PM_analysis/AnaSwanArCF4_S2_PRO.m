%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Script for reading data from Tektronix oscilloscope or CAEN DAQ %%%%%
%%%%%          writing to ASCII file and analyse S1 and S2 signals    %%%%%
%%%%% The final output are the calibrated wvfs and their rms (in pe)  %%%%%
%%%%% for further analysis, together with the main results from the   %%%%%
%%%%% analysis                                                        %%%%%
%%%%%                                                                 %%%%%  
%%%%%                    (DGD 10/Dec/2020)                            %%%%%
%%%%%                  (upgraded 13/Apr/2021)                         %%%%%
%%%%%                                                                 %%%%%
%%%%% NOTES ON Ar-short analysis (S2)                                 %%%%%
%%%%%                                                                 %%%%%
%%%%% Repetition of Xe measurements (see AnaSwanXeShortS1ini)         %%%%%
%%%%% but with Ar                                                     %%%%%
%%%%%                                                                 %%%%%
%%%%% ANALYSIS NOTES (AND SOME TODO's)                                %%%%%
%%%%%                                                                 %%%%%
%%%%% * IMPLEMENT LOOP TO UNPACK ALL DATA                             %%%%%
%%%%% * IMPLEMENT ANALYSIS TO SPLIT S2 and S1 in TWO SAVING TIME FOR  %%%%%
%%%%% LATER ANALYSIS ON DRIFT VELOCITY                                %%%%%
%%%%%                                                                 %%%%%
%%%%% * NICE BECAUSE PHOTON FEEDBACK SEEMS TO APPEAR AS A DECAY CTANT %%%%%
%%%%%                                                                 %%%%%
%%%%% * ASSUME SAME QE for S1 and S2 (might be changed later if       %%%%%
%%%%% needed, based on Florian's S2 data)                             %%%%%
%%%%% * WVF shifted by time at 20% of signal max (configure?)         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;

addpath 'C:\Users\diego\RareEventsGroup Dropbox\HOME_RareEventsGroup\DiegoR\BasicScripts'
addpath 'C:\Users\diego\RareEventsGroup Dropbox\HOME_RareEventsGroup\DiegoR\BasicScripts\_COMMON_SCRIPTS'  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% I. INITIALIZATION (adjust your desired values here)                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic;

NevtsMax = 0;                                                               %Maximum number of events to analyze. Set to zero if not to be used.
DataType = 'CAEN';

% DEFAULT ANALYSIS CUTS                                                     %Might be overriden by Ana_ini
isXcut   = 1; isYcut  = 1;
isA1cut  = 1; isA4cut = 1; 
isT20cut = 1;

AnaSwanArCF4_S2ini;                                                     %Load input file and analysis parameters

isSaveToFile = 1*0; 
DIR_OUTPUT = DIR;                                                           %Only if something is saved

% CUTS FOR S1 analysis
Tmin1_S1 = 0; Tmax1_S1 = 200; 
Tmin2_S1 = 0; Tmax2_S1 = 200; 
Tmin3_S1 = 0; Tmax3_S1 = 200; 
Tmin4_S1 = 0; Tmax4_S1 = 200; 

% SINGLE-PHOTON CALIBRATION                                                 %Do SP analysis with the tails (does not work well for xenon), generally unused... but interesting bckup
isSPanalysis = 0;                                                           
DeltaT_PM1 = 20;  DeltaT_PM2 = 20;  DeltaT_PM3 = 20;  DeltaT_PM4 = 20;      %Time window
Toff_PM1   = 160; Toff_PM2   = 15;  Toff_PM3   = 20;  Toff_PM4   = 160;     %Delay relative to trigger 
sumA_to_ph = 5;                                                             %NOTE: Just approximate, to see the SP-form (6/1/2021). Work with charge!                   

% ANALYSIS PARAMETERS
Rin        = 50;                        %channel input impedance
pedEntries = 200;                       %Time-entries considered in pre-trigger for pedestal computation
isTriggerSetAutomatically = 1;          %Default is one  (no need to specify window            
isSmearedA = 1;                         %smears amplitude
isSmearedT = 1;                         %smears time (FIXME: unused)
NPMTs      = 4;
nSigmaTh   = 3;                         %Number of pedestal sigmas to set threshold in wvf processing (default: 3)

% AUXILLIARY PARAMETERS (for representation, etc)
Arange     = -10:1:4000;                %Range for Amplitude representation [mV]
Qrange     = 0:0.5:4000;                %Range for Charge representation [nC]
Trange     = 0:0.01:5000;               %Time range for fit spectra                        

if(isTriggerSetAutomatically), disp('note: trigger window will be set automatically!'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% II. READOUT (no need to touch))                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(strcmp(DataType,'OSCI')), ReadDataOsci;     end
if(strcmp(DataType,'CAEN')), ReadDataCAEN_PRO; end

if(NevtsMax)
    Nevts = NevtsMax;
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
for k=1:NPMTs
    dataMean(k,:) = mean(squeeze(data_(k,:,:))');                         %#ok<UDIM,SAGROW>
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

%NOTE: threshold set at about x3 the average sigma (this should be tunable)

Amp       = zeros(NPMTs, Nevts); charge    = zeros(NPMTs, Nevts); meanT = zeros(NPMTs, Nevts); 
trise_old = zeros(NPMTs, Nevts); tMin      = zeros(NPMTs, Nevts); trise = zeros(NPMTs, Nevts);
rmsOUT    = zeros(NPMTs, Nevts); stdT      = zeros(NPMTs, Nevts);

for k=1:NPMTs
    Vth(k) = nSigmaTh * mean(pedSTD(k,:));
    
    for i=1:Nevts
        F = data_(k,:,i);  t = time;
        [charge(k,i), Amp(k,i), tMin(k,i), ~, ~, meanT(k,i), trise_old(k,i), ~, rmsOUT(k,i), trise(k,i), ~, ~, ~, ~] = PSA_Nausicaa0(F, t, Vth(k), trigTime);
        
        if(isSmearedA), Amp(k,i)    = Amp(k,i) - ABin/2 + ABin*rand; end
        
        charge(k,i)    = charge(k,i) * TBin / Rin;                   %[mV x ns = pC]
    end
    disp(['PM ', num2str(k), ' digitized']);
end

%Define PM variables as arrays for easier handling

A1      = Amp(1,:);       A2      = Amp(2,:);     A3      = Amp(3,:);     A4      = Amp(4,:);
Q1      = charge(1,:);    Q2      = charge(2,:);  Q3      = charge(3,:);  Q4      = charge(4,:);
Q1norm  = Q1*mean(isPMforNorm)/mean(Q1);  Q2norm  = Q2*mean(isPMforNorm)/mean(Q2);
Q3norm  = Q3*mean(isPMforNorm)/mean(Q3);  Q4norm  = Q4*mean(isPMforNorm)/mean(Q4);  

T1      = meanT(1,:);     T2      = meanT(2,:);   T3      = meanT(3,:);   T4      = meanT(4,:);
stdP1   = pedSTD(1,:);    stdP2   = pedSTD(2,:);  stdP3   = pedSTD(3,:);  stdP4   = pedSTD(4,:);
Base1   = rmsOUT(1,:);    Base2   = rmsOUT(2,:);  Base3   = rmsOUT(3,:);  Base4   = rmsOUT(4,:);

%Trick to get risetime at 20% level.
T20_1   = trise_old(1,:) + tMin(1,:) - trise(1,:);
T20_2   = trise_old(2,:) + tMin(2,:) - trise(2,:);
T20_3   = trise_old(3,:) + tMin(3,:) - trise(3,:);
T20_4   = trise_old(4,:) + tMin(4,:) - trise(4,:);

Yevt  = zeros(1,Nevts); Xevt  = zeros(1,Nevts); Revt  = zeros(1,Nevts);

%Define 'event variables' (normally used for analysis)

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
xaxis(-1, 1); yaxis(-1, 1);
subplot(4,1,3)
hist1D(Xevt, -1:0.01:1);
xlabel('x_{evt} (3-2) [mm]');
subplot(4,1,4)
hist1D(Yevt, -1:0.01:1);
xlabel('y_{evt} (4-1) [mm]');

%PLOT PM variables (no selection)
figure; 
subplot(2,2,1);
plot(Q1, Q1, '.'); title('PM1 vs PM1');
xlabel('Q_{1} [pC]'); ylabel('Q_{1}[pC]'); 
subplot(2,2,2);
plot(Q1, Q2, '.'); title('PM1 vs PM2');
xlabel('Q_{1} [pC]'); ylabel('Q_{2}[pC]');
subplot(2,2,3);
plot(Q1, Q3, '.'); title('PM1 vs PM3');
xlabel('Q_{1} [pC]'); ylabel('Q_{3}[pC]');
subplot(2,2,4);
plot(Q1, Q4, '.'); title('PM1 vs PM4');
xlabel('Q_{1} [pC]'); ylabel('Q_{4}[pC]');

figure;
subplot(2,1,1);
hold on;

for k=1:NPMTs
    if(k==1), hist1D(A1, Arange); lastline('color', 'k');end
    if(k==2), hist1D(A2, Arange); lastline('color', 'r');end
    if(k==3), hist1D(A3, Arange); lastline('color', 'g');end
    if(k==4), hist1D(A4, Arange); lastline('color', 'b');end
end
xaxis(Arange(1),Arange(length(Arange)));
legend('PMT1', 'PMT2', 'PMT3', 'PMT4');
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
N1 = max(hist1D(stdP1, -0.1:0.05:5)); N2 = max(hist1D(stdP2, -0.1:0.05:5));
N3 = max(hist1D(stdP3, -0.1:0.05:5)); N4 = max(hist1D(stdP4, -0.1:0.05:5));

for k=1:NPMTs
    if(k==1), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'k'); lastline('LineStyle', '--'); end
    if(k==2), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'r'); lastline('LineStyle', '--'); end
    if(k==3), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'g'); lastline('LineStyle', '--'); end
    if(k==4), line([Vth(k), Vth(k)], [0,max([N1,N2,N3,N4])]); lastline('color', 'b'); lastline('LineStyle', '--'); end
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4', 'Vth1', 'Vth2', 'Vth3', 'Vth4');
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

Icut   = A1>0;                                                              %unbiased
if(isA1cut),  A1_LOW = 10*mean(stdP1); Icut = Icut & A1>A1_LOW;    end
if(isA4cut),  A4_LOW = 10*mean(stdP4); Icut = Icut & A4>A4_LOW;    end

Icut = Icut & Q1>Q1cutLow & Q1<Q1cutUp;
Icut = Icut & Q2>Q2cutLow & Q2<Q2cutUp;
Icut = Icut & Q3>Q3cutLow & Q3<Q3cutUp;
Icut = Icut & Q4>Q4cutLow & Q4<Q4cutUp;

if(isT20cut),  Icut      = Icut & T20_1<T20_UP & T20_1>T20_LOW; end
if(isXcut)
    rmsX = std(Xevt(Icut)); meanX = mean(Xevt(Icut));
    Xmin = meanX - 1.5*rmsX; Xmax = meanX + 1.5*rmsX;
    Icut = Icut & (Xevt>Xmin) & (Xevt<Xmax); 
end
if(isYcut)
    rmsY = std(Yevt(Icut)); meanY = mean(Yevt(Icut));
    Ymin = meanY - 1.5*rmsY; Ymax = meanY + 1.5*rmsY;
    Icut = Icut & (Yevt>Ymin) & (Yevt<Ymax);  
end

%PLOT PM variables (above threshold)
%PLOT X,Y,R distribution
figure;
subplot(4,1,1)
plot(Xevt(Icut),Yevt(Icut),'.'); lastline('MarkerSize', 1);
if(isXcut)
    line([Xmin, Xmin], [-1,1]); lastline('color', 'r'); lastline('LineStyle', '--');
    line([Xmax, Xmax], [-1,1]); lastline('color', 'r'); lastline('LineStyle', '--');
end
if(isYcut)
    line([-1,1], [Ymin, Ymin]); lastline('color', 'r'); lastline('LineStyle', '--');
    line([-1,1], [Ymax, Ymax]); lastline('color', 'r'); lastline('LineStyle', '--');
end
xlabel('x_{evt} [mm]'); ylabel('y_{evt}[mm]'); title('position after cuts [cut]');
xaxis(-1, 1); yaxis(-1, 1);
subplot(4,1,2)
hist1D(Revt(Icut), 0:0.01:1);
xlabel('R_{evt} [mm]');
subplot(4,1,3)
hist1D(Xevt(Icut), -1:0.01:1);
if(isXcut)
    line([Xmin, Xmin], [0, max(hist1D(Xevt(Icut), -1:0.01:1))]); lastline('color','r');  lastline('LineStyle','--');
    line([Xmax, Xmax], [0, max(hist1D(Xevt(Icut), -1:0.01:1))]); lastline('color','r');  lastline('LineStyle','--');
end
xlabel('X_{evt} [mm]');
subplot(4,1,4)
hist1D(Yevt(Icut), -1:0.01:1);
if(isYcut)
    line([Ymin, Ymin], [0, max(hist1D(Yevt(Icut), -1:0.01:1))]); lastline('color','r');  lastline('LineStyle','--');
    line([Ymax, Ymax], [0, max(hist1D(Yevt(Icut), -1:0.01:1))]); lastline('color','r');  lastline('LineStyle','--');
end
xlabel('Y_{evt} [mm]');

figure;
subplot(2,1,1);
hold on;

for k=1:NPMTs
    if(k==1), hist1D(A1(Icut), Arange); lastline('color', 'k');end
    if(k==2), hist1D(A2(Icut), Arange); lastline('color', 'r');end
    if(k==3), hist1D(A3(Icut), Arange); lastline('color', 'g');end
    if(k==4), hist1D(A4(Icut), Arange); lastline('color', 'b');end
end
if(isA1cut)
    line([A1_LOW, A1_LOW], [0, max(hist1D(A1(Icut), Arange))]); lastline('color','k');  lastline('LineStyle','--');
end
if(isA4cut)
    line([A4_LOW, A4_LOW], [0, max(hist1D(A4(Icut), Arange))]); lastline('color','b');  lastline('LineStyle','--');
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4', 'A1cut', 'A4cut');
xlabel('A[mV]'); ylabel('entries'); title ('Amplitude [cut]'); box on;

subplot(2,1,2);
hold on;
for k=1:NPMTs
    if(k==1), hist1D(Q1(Icut), Qrange); lastline('color', 'k');end
    if(k==2), hist1D(Q2(Icut), Qrange); lastline('color', 'r');end
    if(k==3), hist1D(Q3(Icut), Qrange); lastline('color', 'g');end
    if(k==4), hist1D(Q4(Icut), Qrange); lastline('color', 'b');end
end

line([Q1cutLow, Q1cutLow], [0, max(hist1D(Q1(Icut), Qrange))]); lastline('color','k');  lastline('LineStyle','--');
line([Q2cutLow, Q2cutLow], [0, max(hist1D(Q2(Icut), Qrange))]); lastline('color','r');  lastline('LineStyle','--');
line([Q3cutLow, Q3cutLow], [0, max(hist1D(Q3(Icut), Qrange))]); lastline('color','g');  lastline('LineStyle','--');
line([Q4cutLow, Q4cutLow], [0, max(hist1D(Q4(Icut), Qrange))]); lastline('color','b');  lastline('LineStyle','--');
line([Q1cutUp,  Q1cutUp],  [0, max(hist1D(Q1(Icut), Qrange))]); lastline('color','k');  lastline('LineStyle','--');
line([Q2cutUp,  Q2cutUp],  [0, max(hist1D(Q2(Icut), Qrange))]); lastline('color','r');  lastline('LineStyle','--');
line([Q3cutUp,  Q3cutUp],  [0, max(hist1D(Q3(Icut), Qrange))]); lastline('color','g');  lastline('LineStyle','--');
line([Q4cutUp,  Q4cutUp],  [0, max(hist1D(Q4(Icut), Qrange))]); lastline('color','b');  lastline('LineStyle','--');
   
legend('PMT1', 'PMT2', 'PMT3', 'PMT4', 'Qth1', 'Qth2', 'Qth3', 'Qth4');
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
N1 = max(hist1D(T20_1(Icut), trigWindow(1):trigWindow(2))); N2 = max(hist1D(T20_2(Icut), trigWindow(1):trigWindow(2)));
N3 = max(hist1D(T20_3(Icut), trigWindow(1):trigWindow(2))); N4 = max(hist1D(T20_4(Icut), trigWindow(1):trigWindow(2)));

for k=1:NPMTs
    if(k==1), hist1D(T20_1(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'k'); end
    if(k==2), hist1D(T20_2(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'r'); end
    if(k==3), hist1D(T20_3(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'g'); end
    if(k==4), hist1D(T20_4(Icut), trigWindow(1):trigWindow(2)); lastline('color', 'b'); end
end
if(isT20cut)
    line([T20_LOW, T20_LOW], [0, max([N1,N2,N3,N4])]); lastline('color','k');  lastline('LineStyle','--');
    line([T20_UP,  T20_UP],  [0, max([N1,N2,N3,N4])]); lastline('color','k');  lastline('LineStyle','--');
end
legend('PMT1', 'PMT2', 'PMT3', 'PMT4', 't20cut(PM1)');
xlabel('t[ns]'); ylabel('entries'); title ('time at 20% [cut]'); box on;

%CALCULATE & PLOT (MWVF)
PMT_MWVF_cut  = zeros(NPMTs, dataSize(2));
PMT_MWVF      = zeros(NPMTs, dataSize(2));
MWPC_MWVF_cut = zeros(1, dataSize(2));
MWPC_MWVF     = zeros(1, dataSize(2));

for k=1:(NPMTs+1)
    for i=1:Nevts
        if(Icut(i))
            if(k<5),      PMT_MWVF_cut(k,:)  = PMT_MWVF_cut(k,:)  + data_(k,:,i);
            elseif(k==5 && isQreadout), MWPC_MWVF_cut(1,:) = MWPC_MWVF_cut(1,:) + data(k,:,i);    %MWPC not pedestal corrected
            end
        end
        if(k<5),          PMT_MWVF(k,:)      = PMT_MWVF(k,:)  + data_(k,:,i);
        elseif(k==5 && isQreadout),     MWPC_MWVF(1,:)     = MWPC_MWVF(1,:) + data(k,:,i);        %MWPC not pedestal corrected
        end
    end
end

PMT_MWVF_cut  = PMT_MWVF_cut/length(find(Icut));
PMT_MWVF      = PMT_MWVF/length(Icut);
MWPC_MWVF_cut = MWPC_MWVF_cut/length(find(Icut));
MWPC_MWVF     = MWPC_MWVF/length(Icut);

figure; 
subplot(2,2,1); hold on;
plot(time, PMT_MWVF(1,:));     lastline('color','k');   lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(1,:)); lastline('color','k');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs Ncut)'); box on;

subplot(2,2,2); hold on;
plot(time, PMT_MWVF(2,:));     lastline('color','b');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(2,:)); lastline('color','b');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs Ncut)'); box on;

subplot(2,2,3); hold on;
plot(time, PMT_MWVF(3,:));     lastline('color','r');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(3,:)); lastline('color','r');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs Ncut)'); box on;

subplot(2,2,4); hold on;
plot(time, PMT_MWVF(4,:));     lastline('color','g');  lastline('LineStyle','--');
plot(time, PMT_MWVF_cut(4,:)); lastline('color','g');
logy;
legend('no cuts', 'cuts'); box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF (cut vs Ncut)'); box on;

%CALCULATE & PLOT shifted waveform                                         
PMT_MWVF_cut_shifted = zeros(NPMTs, dataSize(2));                           %Mean
PMT_SWVF_cut_shifted = zeros(NPMTs, dataSize(2));                           %Error of the mean

if(Inew>pedEntries), Inew = pedEntries; end                                 %Do not make larger than pedestal (NOTE: really needed??)
for k=1:NPMTs
    for i=1:Nevts
        if(Icut(i))            
            tf = time(Inew); ti = T20_1(i);
            wvf_i = squeeze(data_(k,:,i));
            wvf_f = TimeAlign(time, wvf_i, ti, tf);
            PMT_MWVF_cut_shifted(k,:) = PMT_MWVF_cut_shifted(k,:) + wvf_f;
            PMT_SWVF_cut_shifted(k,:) = PMT_SWVF_cut_shifted(k,:) + wvf_f.^2;
        end
    end
end

for k=1:NPMTs
    PMT_MWVF_cut_shifted(k,:) = PMT_MWVF_cut_shifted(k,:)/length(find(Icut));
    PMT_SWVF_cut_shifted(k,:) = sqrt(PMT_SWVF_cut_shifted(k,:)/length(find(Icut)) - PMT_MWVF_cut_shifted(k,:).*PMT_MWVF_cut_shifted(k,:));
    PMT_SWVF_cut_shifted(k,:) = PMT_SWVF_cut_shifted(k,:)/sqrt(length(find(Icut)));
end

%Align average waveforms to PM1 maximum (approximate). NOTE: PM1 by default, does not matter except numerically perhaps
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

%Align average waveforms to leave only Inew bins from 0 to max
for k=1:NPMTs    
    I = find(PMT_MWVF_cut(k,:)==max(PMT_MWVF_cut(k,:)));
    tf = time(Inew); ti = time(I);
    PMT_MWVF_cut(k, :) = TimeAlign(time, PMT_MWVF_cut(k, :), ti, tf);
    
    I = find(PMT_MWVF_cut_shifted(k,:)==max(PMT_MWVF_cut_shifted(k,:)));
    tf = time(Inew); ti = time(I);
    PMT_MWVF_cut_shifted(k, :) = TimeAlign(time, PMT_MWVF_cut_shifted(k, :), ti, tf);
    wvf = TimeAlign(time, PMT_SWVF_cut_shifted(k, :), ti, tf);
    wvf(find(wvf==0)) = max(PMT_SWVF_cut_shifted(k, :));                   
    PMT_SWVF_cut_shifted(k, :) = wvf;    
end

%Select regions with S/N larger than 10. [Remove only continuous waveform
%patches left and right from the maximum.]. MOVE TO FUNCTION

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
plot(tim_MWVF_cut1,         PMT_MWVF_cut1);         lastline('color','k');  lastline('LineStyle','--');
legend('time-aligned evt by evt', 'no time-aligned'); logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF'); box on;

subplot(2,2,2); hold on;
plot(tim_MWVF_cut_shifted2, PMT_MWVF_cut_shifted2); lastline('color','b');
plot(tim_MWVF_cut2,         PMT_MWVF_cut2);         lastline('color','b');  lastline('LineStyle','--');
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF'); box on

subplot(2,2,3); hold on;
plot(tim_MWVF_cut_shifted3, PMT_MWVF_cut_shifted3); lastline('color','r');
plot(tim_MWVF_cut3,         PMT_MWVF_cut3);         lastline('color','r');  lastline('LineStyle','--');
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF'); box on;

subplot(2,2,4); hold on;
plot(tim_MWVF_cut_shifted4, PMT_MWVF_cut_shifted4); lastline('color','g');
plot(tim_MWVF_cut4,         PMT_MWVF_cut4);         lastline('color','g');  lastline('LineStyle','--');
logy; box;
xlabel('time [ns]'); ylabel('amplitude [mV]'); title ('mean WVF'); box on;

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
logy; yaxis (0.0001,1.5);
legend('PM1', 'PM2', 'PM3', 'PM4');
xlabel('time [ns]'); ylabel('normalized amplitude '); title ('mean WVF normalized to peak (cut & shifted) '); box on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% VI. SINGLE-PHOTON ANALYSIS (SELDOM USED - BCKUP) %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isSPanalysis)
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
%%%%%%%%%%%%%%%%  VII. OBTAIN CALIBRATED DISTRIBUTIONS     %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wvf1_  = PMT_MWVF_cut_shifted(1,:)*TBin/Rin/(QPM1toPh*TBin);                   %FIXME: This waveform, as given here, is not S/N-suppressed. Better naming needed
wvf2_  = PMT_MWVF_cut_shifted(2,:)*TBin/Rin/(QPM2toPh*TBin);
wvf3_  = PMT_MWVF_cut_shifted(3,:)*TBin/Rin/(QPM3toPh*TBin);
wvf4_  = PMT_MWVF_cut_shifted(4,:)*TBin/Rin/(QPM4toPh*TBin);

swvf1_ = PMT_SWVF_cut_shifted(1,:)*TBin/Rin/(QPM1toPh*TBin);
swvf2_ = PMT_SWVF_cut_shifted(2,:)*TBin/Rin/(QPM2toPh*TBin);
swvf3_ = PMT_SWVF_cut_shifted(3,:)*TBin/Rin/(QPM3toPh*TBin);
swvf4_ = PMT_SWVF_cut_shifted(4,:)*TBin/Rin/(QPM4toPh*TBin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VII.a) S1 WVFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

Ndata      = wvf1_;
sigmaNdata = swvf1_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin1_S1 & Tdata<Tmax1_S1);
sigmaNdata = sigmaNdata(Tdata>Tmin1_S1 & Tdata<Tmax1_S1);
Tdata      = Tdata(Tdata>Tmin1_S1 & Tdata<Tmax1_S1);                                                                                  

if(FitPM1with2Exp && isPM1_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;        %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM1_S1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM1_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM1_S1, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM1with2Exp && isPM1_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM1_S1 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM1_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM1_S1, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,1); hold on;
plot(time, wvf1_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM1with2Exp && isPM1_S1_Fit),      y=spectrum_Exp2Conv(parT_PM1_S1, Trange, 'PM1 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM1with2Exp && isPM1_S1_Fit),   y=spectrum_Exp1Conv(parT_PM1_S1, Trange, 'PM1 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S1');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin1_S1,Tmax1_S1); box;

Q_PM1_Tcut_S1_pe_ = sum(wvf1_(time>Tmin1_S1 & time<Tmax1_S1))*TBin;


%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf2_;
sigmaNdata = swvf2_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin2_S1 & Tdata<Tmax2_S1);
sigmaNdata = sigmaNdata(Tdata>Tmin2_S1 & Tdata<Tmax2_S1);
Tdata      = Tdata(Tdata>Tmin2_S1 & Tdata<Tmax2_S1);                                                                                  

if(FitPM2with2Exp && isPM2_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;        %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM2_S1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM2_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM2_S1, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM2with2Exp && isPM2_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM2_S1 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM2_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM2_S1, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,2); hold on;
plot(time, wvf2_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM2with2Exp && isPM2_S1_Fit),      y=spectrum_Exp2Conv(parT_PM2_S1, Trange, 'PM2 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM2with2Exp && isPM2_S1_Fit),   y=spectrum_Exp1Conv(parT_PM2_S1, Trange, 'PM2 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S1');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin2_S1,Tmax2_S1); box;

Q_PM2_Tcut_S1_pe_ = sum(wvf2_(time>Tmin2_S1 & time<Tmax2_S1))*TBin;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf3_;
sigmaNdata = swvf3_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin3_S1 & Tdata<Tmax3_S1);
sigmaNdata = sigmaNdata(Tdata>Tmin3_S1 & Tdata<Tmax3_S1);
Tdata      = Tdata(Tdata>Tmin3_S1 & Tdata<Tmax3_S1);                                                                                  

if(FitPM3with2Exp && isPM3_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;        %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3_S1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM3_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM3_S1, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM3with2Exp && isPM3_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3_S1 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM3_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM3_S1, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,3); hold on;
plot(time, wvf3_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM3with2Exp && isPM3_S1_Fit),      y=spectrum_Exp2Conv(parT_PM3_S1, Trange, 'PM3 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM3with2Exp && isPM3_S1_Fit),   y=spectrum_Exp1Conv(parT_PM3_S1, Trange, 'PM3 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S1');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin3_S1,Tmax3_S1); box;

Q_PM3_Tcut_S1_pe_ = sum(wvf3_(time>Tmin3_S1 & time<Tmax3_S1))*TBin;           

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf4_;
sigmaNdata = swvf4_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin4_S1 & Tdata<Tmax4_S1);
sigmaNdata = sigmaNdata(Tdata>Tmin4_S1 & Tdata<Tmax4_S1);
Tdata      = Tdata(Tdata>Tmin4_S1 & Tdata<Tmax4_S1);                                                                                  

if(FitPM4with2Exp && isPM4_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;        %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM4_S1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM4_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM4_S1, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM4with2Exp && isPM4_S1_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM4_S1 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM4_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM4_S1, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,4); hold on;
plot(time, wvf4_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM4with2Exp && isPM4_S1_Fit),      y=spectrum_Exp2Conv(parT_PM4_S1, Trange, 'PM4 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM4with2Exp && isPM4_S1_Fit),   y=spectrum_Exp1Conv(parT_PM4_S1, Trange, 'PM4 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S1');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin4_S1,Tmax4_S1); box;

Q_PM4_all_pe_     = sum(wvf4_)*TBin;
Q_PM4_Tcut_S1_pe_ = sum(wvf4_(time>Tmin4_S1 & time<Tmax4_S1))*TBin;           

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% VII.b) S2 WVFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;

Ndata      = wvf1_;
sigmaNdata = swvf1_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin1_S2 & Tdata<Tmax1_S2);
sigmaNdata = sigmaNdata(Tdata>Tmin1_S2 & Tdata<Tmax1_S2);
Tdata      = Tdata(Tdata>Tmin1_S2 & Tdata<Tmax1_S2);                                                                                  

if(FitPM1with2Exp && isPM1_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;       %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM1_S2 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM1_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM1_S2, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM1with2Exp && isPM1_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM1_S2 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM1_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM1_S2, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,1); hold on;
plot(time, wvf1_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM1with2Exp && isPM1_S2_Fit),      y=spectrum_Exp2Conv(parT_PM1_S2, Trange, 'PM1 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM1with2Exp && isPM1_S2_Fit),   y=spectrum_Exp1Conv(parT_PM1_S2, Trange, 'PM1 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S2');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin1_S2,Tmax1_S2);logy; logx; box;

Q_PM1_Tcut_S2_pe_ = sum(wvf1_(time>Tmin1_S2 & time<Tmax1_S2))*TBin;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf2_;
sigmaNdata = swvf2_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin2_S2 & Tdata<Tmax2_S2);
sigmaNdata = sigmaNdata(Tdata>Tmin2_S2 & Tdata<Tmax2_S2);
Tdata      = Tdata(Tdata>Tmin2_S2 & Tdata<Tmax2_S2);                                                                                  

if(FitPM2with2Exp && isPM2_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;       %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM2_S2 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM2_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM2_S2, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM2with2Exp && isPM2_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;       %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM2_S2 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM2_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM2_S2, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,2); hold on;
plot(time, wvf2_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM2with2Exp && isPM2_S2_Fit),      y=spectrum_Exp2Conv(parT_PM2_S2, Trange, 'PM2 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM2with2Exp && isPM2_S2_Fit),   y=spectrum_Exp1Conv(parT_PM2_S2, Trange, 'PM2 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S2');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin2_S2,Tmax2_S2);logy; logx; box;

Q_PM2_Tcut_S2_pe_ = sum(wvf2_(time>Tmin2_S2 & time<Tmax2_S2))*TBin;        

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf3_;
sigmaNdata = swvf3_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin3_S2 & Tdata<Tmax3_S2);
sigmaNdata = sigmaNdata(Tdata>Tmin3_S2 & Tdata<Tmax3_S2);
Tdata      = Tdata(Tdata>Tmin3_S2 & Tdata<Tmax3_S2);                                                                                  

if(FitPM3with2Exp && isPM3_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;        %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3_S2 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM3_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM3_S2, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM3with2Exp && isPM3_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;         %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM3_S2 = lsqnonlin(@fitter_Exp1Conv,par0,lb,ub,options);
    Q_PM3_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM3_S2, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,3); hold on;
plot(time, wvf3_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM3with2Exp && isPM3_S2_Fit),      y=spectrum_Exp2Conv(parT_PM3_S2, Trange, 'PM3 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM3with2Exp && isPM3_S2_Fit),   y=spectrum_Exp1Conv(parT_PM3_S2, Trange, 'PM3 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S2');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin3_S2,Tmax3_S2);logy; logx; box;

Q_PM3_Tcut_S2_pe_ = sum(wvf3_(time>Tmin3_S2 & time<Tmax3_S2))*TBin;

%%%%%%%%%%%%%%%%%%%%%%%%%%% PMT4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ndata      = wvf4_;
sigmaNdata = swvf4_;
Tdata      = time;
Ndata      = Ndata(Tdata>Tmin4_S2 & Tdata<Tmax4_S2);
sigmaNdata = sigmaNdata(Tdata>Tmin4_S2 & Tdata<Tmax4_S2);
Tdata      = Tdata(Tdata>Tmin4_S2 & Tdata<Tmax4_S2);                                                                                  

if(FitPM4with2Exp && isPM4_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude of first exponential
    par0(2)   = 10;         %Time constant of first exponential
    par0(3)   = 10;         %Amplitude of second exponential
    par0(4)   = 100;        %Time constant of second exponential
    par0(5)   = 2;          %Gaussian width of PM response function
    par0(6)   = 250;       %Start of exponentials (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 10000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 0.1;   ub(3)     = 2000;
    lb(4)     = 1;     ub(4)     = 200;
    lb(5)     = 1;     ub(5)     = 5;
    lb(6)     = 245;   ub(6)     = 255;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM4_S2 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM4_fit_pe_  = sum(spectrum_Exp2Conv(parT_PM4_S2, Trange))*(Trange(2)-Trange(1));
elseif(~FitPM4with2Exp && isPM4_S2_Fit)
    clear par0;
    par0(1)   = 10;         %Amplitude
    par0(2)   = 8.5;        %Time constant of exponential
    par0(3)   = 2;          %Gaussian width of PM response function
    par0(4)   = 10;       %Start of exponential (time offset)
    
    lb(1)     = 0.1;   ub(1)     = 1000;
    lb(2)     = 1;     ub(2)     = 40;
    lb(3)     = 1;     ub(3)     = 5;
    lb(4)     = 0;     ub(4)     = 50;
    
    options=optimset('TolX',1e-6,'TolFun',1e-6,'MaxFunEvals', 1000*length(par0),'MaxIter',2000);
    parT_PM4_S2 = lsqnonlin(@fitter_Exp2Conv,par0,lb,ub,options);
    Q_PM4_fit_pe_  = sum(spectrum_Exp1Conv(parT_PM4_S2, Trange))*(Trange(2)-Trange(1));
end

subplot(2,2,4); hold on;
plot(time, wvf4_, '+', 'Markersize', 10);  lastline('color', 'k');

if(FitPM4with2Exp && isPM4_S2_Fit),      y=spectrum_Exp2Conv(parT_PM4_S2, Trange, 'PM4 (final)'); xaxis(min(Trange),max(Trange));
elseif(~FitPM4with2Exp && isPM4_S2_Fit),   y=spectrum_Exp1Conv(parT_PM4_S2, Trange, 'PM4 (final)'); xaxis(min(Trange),max(Trange));
end

xlabel('time[ns]'); ylabel('spectrum [pe/ns]'); title('S2');
lastline('LineStyle','--');
yaxis(max(Ndata)/1000,max(Ndata)*1.2); xaxis(Tmin4_S2,Tmax4_S2);logy; logx; box;

Q_PM4_Tcut_S2_pe_ = sum(wvf4_(time>Tmin4_S2 & time<Tmax4_S2))*TBin; 

figure; hold on;

plot(time/1000, (MWPC_MWVF     - MWPC_MWVF(1)),     '-');
plot(time/1000, (MWPC_MWVF_cut - MWPC_MWVF_cut(1)), 'r-');
xlabel('time[us]'); ylabel('A[mV]'); legend('uncut','cut'); title('average MWPC signal'); box;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% VII. CALCULATE W %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ealpha   = 5.5*1e+6;                    %Energy of alpha particle.
Omega_S1 = 0.0085 * 1.72;               %Makes average between filter and no filter (0.00875 and 0.0081) for nominal (Xe-paper) geometry. 
                                        %Introduces a factor 1.72 to reproduce Xe-Wsc in S1-analysis (simulation by Angela gives ~2.7, 
                                        %but not in exact same geometry. It should be less, perhaps 2).
Omega_S2 = 2 * 0.0085 * 1.72;           %Approximate value, assuming r^2 dependence (the increase should be probably less)

QE_PM1   = 0.21;                        %PM7378 efficiency (typical -> take the same as PM2)
if (isAr), QE_PM1   = 0.185; end        %PM7378 efficiency for Ar 3rd continuum

QE_PM2   = 0.21;                        %PM5070 efficiency (from Pablo's spectrum: 1%CF4->0.21, 10%CF4->0.2095). Done for S1.
if (isAr), QE_PM2   = 0.049; end        %PM7378 efficiency for Ar 3rd continuum with 250-400 filter

QE_PM3   = 0.0802;                      %PM5070 efficiency (from Pablo's spectrum: 1%CF4->0.0816, 10%CF4->0.0788). Done for S1.
if (isAr), QE_PM3   = 6.9770e-05; end   %PM5070 efficiency for Ar 3rd continuum with 400-700 filter

Omega_S1_wfilter = 0.0085 * 1.72;       %Assume the same acceptance for filtered and unfiltered
Omega_S2_wfilter = 2 * 0.0085 * 1.72;   %Assume the same acceptance for filtered and unfiltered

QE_400nm_PM4 = 0.25;                    %Quantum efficiency of PM4 at 400nm
TPB_WLE      = 0.2;                     %WLE of TPB
OmegaTPB     = 0.3;                     %Approximate solid angle of TPB
QE_PM4       = TPB_WLE * OmegaTPB * QE_400nm_PM4; %Quantum efficiency of PM4 to UV light

Transp  = 0.936*(0.97^2)*0.99;      %Transparency of anode (thin) and 3 additional planes.
WI      = 26.4;                     %energy to create an electron-ion pair for argon

nph_prod_S1_PM1 = Q_PM1_Tcut_S1_pe_ /(Omega_S1*QE_PM1*Transp);
W_S1_PM1        = Ealpha/nph_prod_S1_PM1;

nph_prod_S2_PM1 = Q_PM1_Tcut_S2_pe_ /(Omega_S2*QE_PM1*Transp);
W_S2_PM1        = Ealpha/nph_prod_S2_PM1;

nph_prod_S1_PM2 = Q_PM2_Tcut_S1_pe_ /(Omega_S1*QE_PM2*Transp);
W_S1_PM2        = Ealpha/nph_prod_S1_PM2;

nph_prod_S2_PM2 = Q_PM2_Tcut_S2_pe_ /(Omega_S2*QE_PM2*Transp);
W_S2_PM2        = Ealpha/nph_prod_S2_PM2;

nph_prod_S1_PM3 = Q_PM3_Tcut_S1_pe_ /(Omega_S1*QE_PM3*Transp);
W_S1_PM3        = Ealpha/nph_prod_S1_PM3;

nph_prod_S2_PM3 = Q_PM3_Tcut_S2_pe_ /(Omega_S2*QE_PM3*Transp);
W_S2_PM3        = Ealpha/nph_prod_S2_PM3;

nph_prod_S1_PM4 = Q_PM4_Tcut_S1_pe_ /(Omega_S1*QE_PM4*Transp);
W_S1_PM4        = Ealpha/nph_prod_S1_PM4;

nph_prod_S2_PM4 = Q_PM4_Tcut_S2_pe_ /(Omega_S2*QE_PM4*Transp);
W_S2_PM4        = Ealpha/nph_prod_S2_PM4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% VII. SUMMARY: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('                                                                   ');
disp('------------------------   summary   ------------------------------');
disp('                                                                   ');

disp('                PM1 R7378 (no filter)                    ');

disp(['W_{s1} [eV]  (250-400nm) =            ', num2str(W_S1_PM2)]);
disp(['W_{s1} [eV]  (400-700nm) =            ', num2str(W_S1_PM3)]);
disp(['optical gain (400-700nm) =            ', num2str(WI/W_S2_PM3)]);
disp(['optical gain (250-400nm) =            ', num2str(WI/W_S2_PM2)]);

fclose('all');
toc;

mosaic;
return;

%L1517, L1463
