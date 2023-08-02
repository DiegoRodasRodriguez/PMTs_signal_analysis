%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% RUN CORRESPONDING TO DATA TAKEN IN MARCH 2021 IN TPC WITH 0.7CM DRIFT   %
% DISTANCE AND MWPC INSTALLED. THERE IS NEITHER TEFLON REFLECTOR NOR FIELD%
% SHAPERS. THERE IS AN ALPHA SOURCE, WITH 3 PMT R7378 and 1 5070 PMT (3)  %
% FACING IT. PM1 HAS NO FILTER, PM2 HAS A 250-400nm FILTER, PM3 HAS A     %
% 250-700nm FILTER AND PM4 HAS TPB COATING ON IT                          %
%                                                                         %
% THERE ARE FILES WITH S1, AND WITH S2, AND S2 FILES CONTAIN S1 TOO       %
% (SO THEY ALLOW TO ESTIMATE THE DRIFT VELOCITY)                          %
%                                                                         %
% THIS INITIALITATION FILE CORRESPONDS TO S1-ANALYSIS ONLY                %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nsave = 12000;
FILE = 'wave';

%Usually don't need to touch, but for very narrow WVFs this cut can lead to
%insufficient points to fit. PROBALY ISSUE A WARNING IN FINAL OUTPUT.
StoNcut_PM1 = 10; isPM1_StoNcut = 1*0;    
StoNcut_PM2 = 10; isPM2_StoNcut = 1*0;
StoNcut_PM3 = 10; isPM3_StoNcut = 1*0;
StoNcut_PM4 = 10; isPM4_StoNcut = 1*0;

isPMforNorm  = 1;                                                            %PM for relative normalization in order to estimate position (1 by default)
isPMforShift = 1;                                                            %Shift all average waveforms to match PM1
TakeAfterPulsingFromExt = 0;

QPM1toPh = 0.88*(1.15/1.35)^7.63;   % 0.88 +- 0.07
QPM2toPh = 0.45*(1.35/1.35)^7.63;   % 0.45 +- 0.03
QPM3toPh = 0.45*(1.25/1.35)^7.63;   % 0.45 +- 0.04
QPM4toPh = 0.42*(1.25/1.35)^7.63;   % 0.42 +- 0.04

isPM3Fit = 1;   
nCh      = 5;                                                                %For acquisition

isAr = 0;
isQreadout = 1;                                                              %For reading out anode

%%%%%%%%%%%% FILES %%%%%%%%%%%%%%%%%%%%%
                               
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10bar campaign (S1+S2) in SWAN %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                           Ar/CF4                                %%%%%  
%%%%% NOTE: afterpulsing corrected from Ar/CF4 at highest fraction    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FitPM1with2Exp = 0; isPM1_S1_Fit = 0; isPM1_S2_Fit = 0;
FitPM2with2Exp = 0; isPM2_S1_Fit = 0; isPM2_S2_Fit = 0;
FitPM3with2Exp = 0; isPM3_S1_Fit = 0; isPM3_S2_Fit = 0;
FitPM4with2Exp = 0; isPM4_S1_Fit = 0; isPM4_S2_Fit = 0;

%%%%%%%%%%%%%%%%%
%%%   0.2%   %%%%
%%%%%%%%%%%%%%%%%
% Inew     = 165;                                            %Bin number where to realign waveform
% Tmin1_S2 = 640; Tmax1_S2 = 750; 
% Tmin2_S2 = 550; Tmax2_S2 = 750; 
% Tmin3_S2 = 550; Tmax3_S2 = 750; 
% Tmin4_S2 = 600; Tmax4_S2 = 750; 

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3000_V\'; fCF4 = 0.2;
%Q1cutLow = 5; Q1cutUp = 30; Q2cutLow = 0; Q2cutUp = 20; Q3cutLow = 11.5; Q3cutUp = 33.5; Q4cutLow = 30; Q4cutUp = 90; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3100_V\'; fCF4 = 0.2; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3200_V\'; fCF4 = 0.2;
%Q1cutLow = 12; Q1cutUp = 40; Q2cutLow = 2; Q2cutUp = 20; Q3cutLow = 18; Q3cutUp = 50; Q4cutLow = 40; Q4cutUp = 110; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3300_V\'; fCF4 = 0.2; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3400_V\'; fCF4 = 0.2;
%Q1cutLow = 20; Q1cutUp = 45; Q2cutLow = 0; Q2cutUp = 20; Q3cutLow = 30; Q3cutUp = 70; Q4cutLow = 50; Q4cutUp = 160; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3500_V\'; fCF4 = 0.2; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3600_V\'; fCF4 = 0.2;
%Q1cutLow = 30; Q1cutUp = 65; Q2cutLow = 0; Q2cutUp = 25; Q3cutLow = 50; Q3cutUp = 100; Q4cutLow = 80; Q4cutUp = 250; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3700_V\'; fCF4 = 0.2; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\29March_99_8_Ar_0_2_CF4\3800_V\'; fCF4 = 0.2; %undershoot in PM1!
%Q1cutLow = 10; Q1cutUp = 30; Q2cutLow = 5; Q2cutUp = 30; Q3cutLow = 70; Q3cutUp = 160; Q4cutLow = 125; Q4cutUp = 350; T20_LOW = 2770; T20_UP = 2870;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   0.6%  (copied from 2%. Are these files actually analyzed?    %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3485_V\'; fCF4 = 0.6;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3585_V\'; fCF4 = 0.6; %TBunpacked
%Q1cutLow = 40; Q1cutUp = 90; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 80; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3685_V\'; fCF4 = 0.6; 
%Q1cutLow = 40; Q1cutUp = 90; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 80; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3785_V\'; fCF4 = 0.6; %TBunpacked
%Q1cutLow = 50; Q1cutUp = 120; Q2cutLow = 5; Q2cutUp = 40; Q3cutLow = 100; Q3cutUp = 250; Q4cutLow = 60; Q4cutUp = 190; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3885_V\'; fCF4 = 0.6; 
%Q1cutLow = 60; Q1cutUp = 150; Q2cutLow = 10; Q2cutUp = 50; Q3cutLow = 150; Q3cutUp = 270; Q4cutLow = 70; Q4cutUp = 250; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3985_V\'; fCF4 = 0.6; %TBunpacked
% Q1cutLow = 90; Q1cutUp = 190; Q2cutLow = 10; Q2cutUp = 60; Q3cutLow = 200; Q3cutUp = 400; Q4cutLow = 100; Q4cutUp = 300; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4085_V\'; fCF4 = 0.6; 
%Q1cutLow = 125; Q1cutUp = 250; Q2cutLow = 15; Q2cutUp = 80; Q3cutLow = 250; Q3cutUp = 500; Q4cutLow = 150; Q4cutUp = 350; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4185_V\'; fCF4 = 0.6; %TBunpacked
% Q1cutLow = 150; Q1cutUp = 350; Q2cutLow = 20; Q2cutUp = 100; Q3cutLow = 200; Q3cutUp = 800; Q4cutLow = 150; Q4cutUp = 600; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4285_V\'; fCF4 = 0.6; 
% Q1cutLow = 125; Q1cutUp = 500; Q2cutLow = 30; Q2cutUp = 120; Q3cutLow = 400; Q3cutUp = 1000; Q4cutLow = 250; Q4cutUp = 700; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4385_V\'; fCF4 = 0.6; %TBunpacked
% Q1cutLow = 125; Q1cutUp = 700; Q2cutLow = 50; Q2cutUp = 170; Q3cutLow = 500; Q3cutUp = 1500; Q4cutLow = 250; Q4cutUp = 1000; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4485_V\'; fCF4 = 0.6; 
%Q1cutLow = 400; Q1cutUp = 1100; Q2cutLow = 50; Q2cutUp = 250; Q3cutLow = 1100; Q3cutUp = 1900; Q4cutLow = 450; Q4cutUp = 1400; T20_LOW = 2770; T20_UP = 2870; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4530_V\'; fCF4 = 0.6; 
%Q1cutLow = 600; Q1cutUp = 1700; Q2cutLow = 70; Q2cutUp = 350; Q3cutLow = 1200; Q3cutUp = 2400; Q4cutLow = 750; Q4cutUp = 2000; T20_LOW = 2770; T20_UP = 2870; 

%%%%%%%%%%%%%%%%%
%%%   1.0%   %%%%
%%%%%%%%%%%%%%%%%
Inew     = 100;                                            %Bin number where to realign waveform (depends on gas // drift velocity)
Tmin1_S2 = 0; Tmax1_S2 = 1000;                             %These numbers will require readjustment
Tmin2_S2 = 0; Tmax2_S2 = 1000; 
Tmin3_S2 = 0; Tmax3_S2 = 1000; 
Tmin4_S2 = 0; Tmax4_S2 = 1000; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3300_V\'; fCF4 = 1;
%Q1cutLow = 20; Q1cutUp = 60; Q2cutLow = 2; Q2cutUp = 25; Q3cutLow = 50; Q3cutUp = 110; Q4cutLow = 30; Q4cutUp = 90; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3400_V\'; fCF4 = 1;
%Q1cutLow = 40; Q1cutUp = 90; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 80; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3500_V\'; fCF4 = 1;
%Q1cutLow = 40; Q1cutUp = 90; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 80; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3600_V\'; fCF4 = 1;
%Q1cutLow = 50; Q1cutUp = 120; Q2cutLow = 5; Q2cutUp = 40; Q3cutLow = 100; Q3cutUp = 250; Q4cutLow = 60; Q4cutUp = 190; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3700_V\'; fCF4 = 1;
%Q1cutLow = 60; Q1cutUp = 150; Q2cutLow = 10; Q2cutUp = 50; Q3cutLow = 150; Q3cutUp = 270; Q4cutLow = 70; Q4cutUp = 250; T20_LOW = 2770; T20_UP = 2870;
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3800_V\'; fCF4 = 1;
% Q1cutLow = 90; Q1cutUp = 190; Q2cutLow = 10; Q2cutUp = 60; Q3cutLow = 200; Q3cutUp = 400; Q4cutLow = 100; Q4cutUp = 300; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\3900_V\'; fCF4 = 1;
%Q1cutLow = 125; Q1cutUp = 250; Q2cutLow = 15; Q2cutUp = 80; Q3cutLow = 250; Q3cutUp = 500; Q4cutLow = 150; Q4cutUp = 350; T20_LOW = 2770; T20_UP = 2870;
DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\4000_V\'; fCF4 = 1; %W_uv = 4257eV; W_vis = 630eV; gain_uv = 0.1; gain_vis = 4.7;
Q1cutLow = 150; Q1cutUp = 350; Q2cutLow = 20; Q2cutUp = 100; Q3cutLow = 200; Q3cutUp = 800; Q4cutLow = 150; Q4cutUp = 600; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\4100_V\'; fCF4 = 1;
%Q1cutLow = 125; Q1cutUp = 500; Q2cutLow = 30; Q2cutUp = 120; Q3cutLow = 400; Q3cutUp = 1000; Q4cutLow = 250; Q4cutUp = 700; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\4200_V\'; fCF4 = 1;
%Q1cutLow = 125; Q1cutUp = 700; Q2cutLow = 50; Q2cutUp = 170; Q3cutLow = 500; Q3cutUp = 1500; Q4cutLow = 250; Q4cutUp = 1000; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\4300_V\'; fCF4 = 1; %No charge amplifier
%Q1cutLow = 400; Q1cutUp = 1100; Q2cutLow = 50; Q2cutUp = 250; Q3cutLow = 1100; Q3cutUp = 1900; Q4cutLow = 450; Q4cutUp = 1400; T20_LOW = 2770; T20_UP = 2870; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_99_Ar_1_CF4\4350_V\'; fCF4 = 1;
%Q1cutLow = 600; Q1cutUp = 1700; Q2cutLow = 70; Q2cutUp = 350; Q3cutLow = 1200; Q3cutUp = 2400; Q4cutLow = 750; Q4cutUp = 2000; T20_LOW = 2770; T20_UP = 2870; 

%%%%%%%%%%%%%%%%%
%%%   2.0%   %%%%
%%%%%%%%%%%%%%%%%

% Inew     = 100;                                            %Bin number where to realign waveform
% Tmin1_S2 = 370; Tmax1_S2 = 550; 
% Tmin2_S2 = 350; Tmax2_S2 = 520; 
% Tmin3_S2 = 350; Tmax3_S2 = 440; 
% Tmin4_S2 = 350; Tmax4_S2 = 550; 

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3485_V\'; fCF4 = 2.0;
%Q1cutLow = 30; Q1cutUp = 65; Q2cutLow = 0; Q2cutUp = 30; Q3cutLow = 30; Q3cutUp = 60; Q4cutLow = 30; Q4cutUp = 90; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3585_V\'; fCF4 = 2.0; %TBunpacked
%Q1cutLow = 40; Q1cutUp = 90; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 80; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3685_V\'; fCF4 = 2.0; 
%Q1cutLow = 40; Q1cutUp = 100; Q2cutLow = 4; Q2cutUp = 35;  Q3cutLow = 40; Q3cutUp = 180; Q4cutLow = 50; Q4cutUp = 130; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3785_V\'; fCF4 = 2.0; %TBunpacked
%Q1cutLow = 50; Q1cutUp = 120; Q2cutLow = 5; Q2cutUp = 40; Q3cutLow = 100; Q3cutUp = 250; Q4cutLow = 60; Q4cutUp = 190; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3885_V\'; fCF4 = 2.0; 
%Q1cutLow = 80; Q1cutUp = 160; Q2cutLow = 10; Q2cutUp = 50; Q3cutLow = 170; Q3cutUp = 320; Q4cutLow = 80; Q4cutUp = 220; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\3985_V\'; fCF4 = 2.0; %TBunpacked
%Q1cutLow = 90; Q1cutUp = 190; Q2cutLow = 10; Q2cutUp = 60; Q3cutLow = 200; Q3cutUp = 400; Q4cutLow = 100; Q4cutUp = 300; T20_LOW = 2770; T20_UP = 2870;
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4085_V\'; fCF4 = 2.0; 
% Q1cutLow = 120; Q1cutUp = 300; Q2cutLow = 15; Q2cutUp = 80; Q3cutLow = 300; Q3cutUp = 550; Q4cutLow = 125; Q4cutUp = 400; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4185_V\'; fCF4 = 2.0; %TBunpacked
% Q1cutLow = 150; Q1cutUp = 350; Q2cutLow = 20; Q2cutUp = 100; Q3cutLow = 200; Q3cutUp = 800; Q4cutLow = 150; Q4cutUp = 600; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4285_V\'; fCF4 = 2.0; 
%Q1cutLow = 200; Q1cutUp = 550; Q2cutLow = 30; Q2cutUp = 150; Q3cutLow = 500; Q3cutUp = 1200; Q4cutLow = 250; Q4cutUp = 800; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4385_V\'; fCF4 = 2.0; %TBunpacked
% Q1cutLow = 125; Q1cutUp = 700; Q2cutLow = 50; Q2cutUp = 170; Q3cutLow = 500; Q3cutUp = 1500; Q4cutLow = 250; Q4cutUp = 1000; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4485_V\'; fCF4 = 2.0; 
%Q1cutLow = 500; Q1cutUp = 1400; Q2cutLow = 60; Q2cutUp = 270; Q3cutLow = 1100; Q3cutUp = 2400; Q4cutLow = 400; Q4cutUp = 1400; T20_LOW = 2770; T20_UP = 2870; 
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\26March_98_Ar_2_CF4\4530_V\'; fCF4 = 2.0; 
% Q1cutLow = 600; Q1cutUp = 1800; Q2cutLow = 100; Q2cutUp = 350; Q3cutLow = 1500; Q3cutUp = 2600; Q4cutLow = 700; Q4cutUp = 2000; T20_LOW = 2770; T20_UP = 2870; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   1.0%CF4 + 0.5%CH4  %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inew     = 100;                                            %Bin number where to realign waveform
% Tmin1_S2 = 370; Tmax1_S2 = 550; 
% Tmin2_S2 = 350; Tmax2_S2 = 520;  
% Tmin3_S2 = 350; Tmax3_S2 = 440; 
% Tmin4_S2 = 350; Tmax4_S2 = 550; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\5April_0_5_CH4\3600_V\'; fCF4 = 1; fCH4 = 0.5;
%Q1cutLow = 100; Q1cutUp = 225; Q2cutLow = 10; Q2cutUp = 70; Q3cutLow = 220; Q3cutUp = 500; Q4cutLow = 90; Q4cutUp = 240; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\5April_0_5_CH4\3700_V\'; fCF4 = 1; fCH4 = 0.5; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\5April_0_5_CH4\3800_V\'; fCF4 = 1; fCH4 = 0.5;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\5April_0_5_CH4\3900_V\'; fCF4 = 1; fCH4 = 0.5; %TBunpacked
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\5April_0_5_CH4\4000_V\'; fCF4 = 1; fCH4 = 0.5;
%Q1cutLow = 150; Q1cutUp = 600; Q2cutLow = 10; Q2cutUp = 140; Q3cutLow = 600; Q3cutUp = 1200; Q4cutLow = 150; Q4cutUp = 600; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6April_0_5_CH4\4100_V\'; fCF4 = 1; fCH4 = 0.5;
%Q1cutLow = 150; Q1cutUp = 600; Q2cutLow = 10; Q2cutUp = 140; Q3cutLow = 600; Q3cutUp = 1200; Q4cutLow = 150; Q4cutUp = 600; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6April_0_5_CH4\4200_V\'; fCF4 = 1; fCH4 = 0.5;
%Q1cutLow = 150; Q1cutUp = 600; Q2cutLow = 25; Q2cutUp = 350; Q3cutLow = 600; Q3cutUp = 1600; Q4cutLow = 350; Q4cutUp = 800; T20_LOW = 2770; T20_UP = 2870;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   2%CH4              %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inew     = 100;                                            %Bin number where to realign waveform
%Tmin1_S2 = 370; Tmax1_S2 = 550; 
%Tmin2_S2 = 350; Tmax2_S2 = 520; 
%Tmin3_S2 = 350; Tmax3_S2 = 440; 
%Tmin4_S2 = 350; Tmax4_S2 = 550; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3300_V\'; fCH4 = 2;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3400_V\'; fCF4 = 2; fCH4 = 2; %T-window crazy
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3500_V\'; fCF4 = 2; fCH4 = 2; %T-window crazy
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3600_V\'; fCF4 = 2; fCH4 = 2;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3700_V\'; fCF4 = 2; fCH4 = 2;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3800_V\'; fCF4 = 2; fCH4 = 2;
%Q1cutLow = 20*0; Q1cutUp = 60*10; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 50*0; Q3cutUp = 110*10; Q4cutLow = 30*0; Q4cutUp = 90*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\6april_2_CH4_98_Ar_NO_CF4\3900_V\'; fCF4 = 2; fCH4 = 2;
%Q1cutLow = 100; Q1cutUp = 900; Q2cutLow = 2*0; Q2cutUp = 25*10; Q3cutLow = 70; Q3cutUp = 200; Q4cutLow = 70; Q4cutUp = 200; T20_LOW = 2770; T20_UP = 2870;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   pure argon         %%%% Here gain higher even than in Ar/CF4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inew     = 100;                                            %Bin number where to realign waveform
%Tmin1_S2 = 0; Tmax1_S2 = 550; 
%Tmin2_S2 = 0; Tmax2_S2 = 520; 
%Tmin3_S2 = 0; Tmax3_S2 = 440; 
%Tmin4_S2 = 0; Tmax4_S2 = 550; 
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\7April_100_Ar\3000_V\'; fCH4 = 0;
%Q1cutLow = 100*0; Q1cutUp = 220*10; Q2cutLow = 10*0; Q2cutUp = 70*10; Q3cutLow = 220*0; Q3cutUp = 500*10; Q4cutLow = 90*0; Q4cutUp = 240*10; T20_LOW = 2770; T20_UP = 2870;
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\7April_100_Ar\3100_V\'; fCH4 = 0;
% Q1cutLow = 100*0; Q1cutUp = 220*10; Q2cutLow = 10*0; Q2cutUp = 70*10; Q3cutLow = 220*0; Q3cutUp = 500*10; Q4cutLow = 90*0; Q4cutUp = 240*10; T20_LOW = 2770; T20_UP = 2870;
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\7April_100_Ar\3200_V\'; fCH4 = 0;   %ANALYZED RUN
% Q1cutLow = 100*0; Q1cutUp = 220*10; Q2cutLow = 10*0; Q2cutUp = 70*10; Q3cutLow = 220*0; Q3cutUp = 500*10; Q4cutLow = 90*0; Q4cutUp = 240*10; T20_LOW = 2770; T20_UP = 2870;
% DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\7April_100_Ar\3300_V\'; fCH4 = 0;
% Q1cutLow = 100*0; Q1cutUp = 220*10; Q2cutLow = 10*0; Q2cutUp = 70*10; Q3cutLow = 220*0; Q3cutUp = 500*10; Q4cutLow = 90*0; Q4cutUp = 240*10; T20_LOW = 2770; T20_UP = 2870;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\7April_100_Ar\3400_V\'; fCH4 = 0;
%Q1cutLow = 100*0; Q1cutUp = 220*10; Q2cutLow = 10*0; Q2cutUp = 70*10; Q3cutLow = 220*0; Q3cutUp = 500*10; Q4cutLow = 90*0; Q4cutUp = 240*10; T20_LOW = 2770; T20_UP = 2870;
%isAr = 1;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%10bar campaign (S1+S2) in BAT  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                           Ar/CF4 (0.7/99.3)                     %%%%%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inew     = 100;                                            %Bin number where to realign waveform (depends on gas // drift velocity)
%Tmin1_S2 = 0; Tmax1_S2 = 1000;                             %These numbers will require readjustment
%Tmin2_S2 = 0; Tmax2_S2 = 1000; 
%Tmin3_S2 = 0; Tmax3_S2 = 1000; 
%Tmin4_S2 = 0; Tmax4_S2 = 1000; 

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\BAT\Data\2023\RUN_Nov2022\21_11_2022\PMTs\';
%Q1cutLow = 0; Q1cutUp = 350; Q2cutLow = 0; Q2cutUp = 100; Q3cutLow = 0; Q3cutUp = 200; Q4cutLow = 0; Q4cutUp = 200; T20_LOW = 3700; T20_UP = 3870; isXcut = 0; isYcut = 0; isQreadout = 0;
