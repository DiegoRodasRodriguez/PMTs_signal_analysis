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

StoNcut_PM1 = 10; isPM1_StoNcut = 1*0;                                      %For suppression of regions with low S/N levels
StoNcut_PM2 = 10; isPM2_StoNcut = 1*0;                                      %For suppression of regions with low S/N levels
StoNcut_PM3 = 10; isPM3_StoNcut = 1*0;                                      %For suppression of regions with low S/N levels
StoNcut_PM4 = 10; isPM4_StoNcut = 1*0;                                      %For suppression of regions with low S/N levels

TakeAfterPulsingFromExt = 0;                                                %Take afterpulsing from existing template
isXe   = 0;                                                                 %Includes QE for Xe emissions
isCorr = 0;                                                                 %Corrects for triplet quenching in case of Xe
isAr   = 0;                                                                 %Includes QE for Ar emissions

%REFERENCE ANALYSIS (CALIBRATION USED FOR XE PAPER: EXPONENTIAL)
%QPM1toPh = 1.02*exp(-7.214*0.2);
%QPM2toPh = 0.51*exp(-7.124*0.0);
%QPM3toPh = 0.97*exp(-6.966*0.1);
%QPM4toPh = 0.50*exp(-7.1*0.1);

%REANALYSIS (IMPROVED CALIBRATION: POWER LAW)
QPM1toPh = 0.88*(1.15/1.35)^7.63;   % 0.88 +- 0.07
QPM2toPh = 0.45*(1.35/1.35)^7.63;   % 0.45 +- 0.03
QPM3toPh = 0.45*(1.25/1.35)^7.63;   % 0.45 +- 0.04
QPM4toPh = 0.42*(1.25/1.35)^7.63;   % 0.42 +- 0.04

%%%%%%%%%%%% FILES %%%%%%%%%%%%%%%%%%%%%
                         
%%%%%%%%%%%%%%%%
%10bar campaign%
%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%
%%%%% xenon %%%%
%%%%%%%%%%%%%%%%

%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\10barNEWPMTS\1000V\';
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\10barNEWPMTS\1000Vafter1h\';
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\13AprilXenon\'; %disastrous purity!
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%TakeAfterPulsingFromExt = 1; DIR_OUTPUTa = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\10barNEWPMTS\';
%isCorr = 1;
%isXe   = 1;

% PM1 (2ND CONTINUUM) -> OK IF TUNING A BIT GEOM EFF (SOME 50%).  TIME CTANT FINE
% PM2 (THIRD CONTINUUM +250-400 filter) -> SEES ABOUT x1/10 of expected. TIME COTANT WRONG -> PROBLEM WITH PM??
% PM3 (THIRD CONTINUUM +400-700 filter) -> 1% OK! TME CTANT NOT TOO BAD.
% PM4 (2ND CONTINUUM +400-700 filter)   -> ALL FINE IF ASSUMING 20% TPB_WLSE

%%%%%%%%%%%%%%%%
%%%%% argon %%%%
%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\1000V\'; %NOTE: PM1 biased in amplitude by ~x1.5 //nBr clearly visible!
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\575V\';  %NOTE: PM1 biased in amplitude by ~x1.2
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\306V\';  %NOTE: PM1 biased in amplitude by ~x1.3
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\0V\';    %NOTE: PM1 biased in amplitude by ~x1.3
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\115V\';  %NOTE: PM1 biased in amplitude by ~x1.3
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\192V\';  %NOTE: PM1 biased in amplitude by ~x1.3
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%TakeAfterPulsingFromExt = 1; DIR_OUTPUTa = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon\NewSetup\10bar\AfterpulsesSaved\';
%Temporary afterpulsing library using PM1 as a reference. Take it in the future from ArCF4 at high quencher rates.
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%isA1cut = 1; isA4cut = 1;
%isAr   = 1;

%NOTES (with PM1 biased in amplitude by ~x2):
% PM1 (THIRD CONTINUUM) -> SEES JUST 7% (2.5 times less than SANTORELLI'S PAPER BUT threshold-biased). TIME CTANT FINE.
% PM2 (THIRD CONTINUUM +250-400 filter) -> SEES ABOUT 1/2 of PM1 (4% when extrapolated to entire area). TIME CTANT FINE
%IN FACT, IF REMOVING PM1 BIAS THE AGREEMENT BETWEEN PM1 AND PM2 AT AROUND 3-5% SEEMS MUCH BETTER (x4-x6 DISAGREEMENT WITH SANTORELLI!)
% PM3 (THIRD CONTINUUM +400-700 filter) -> TOO LITTLE COVERAGE. OVERCORRECTS!
% SUGGESTS THAT PM2 PROBABLY NOT TOO WRONG ALTHOUGH OVERALL x3-x6 LESS
% LIGHT THAN SANTORELLI. SUGGESTS THAT Ar-states easier to quench than Xe-ones

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%                           Ar/CF4                                %%%%%  
%%%%% NOTE: afterpulsing corrected from Ar/CF4 at highest fraction    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
%%%    0.1%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\576V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\756V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_9_Ar_0_1_CF4\756Vfinal\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;

%%%%%%%%%%%%%%%%%
%%%    0.2%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%
% PM2 and PM3 have been exchanged during the readout! Patched by simply swapping the channels in main

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\576V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_8_Ar_0_2_CF4\756V\';

%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;

%%%%%%%%%%%%%%%%%
%%%    0.5%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\576V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\762V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_5_Ar_0_5_CF4\762Vfinal\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;

%%%%%%%%%%%%%%%%%
%%%    1.0%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\306V\'; %PM1 amplitude biased
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\576V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\756V\';
DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\756Vfinal\';
useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
TakeAfterPulsingFromExt = 1; DIR_OUTPUTa = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\';    
isA1cut = 1; isA4cut = 1;
QPM2toPh = 0.45*(1.25/1.35)^7.63;   % 0.45 +- 0.03                          %NOTE: correct for small mistake on PM HV during run

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4/NewSetup/10bar/18March_99_Ar_1_CF4/primary_230V/';
%useAfterPulsingTemplate=1*0; storeAfterPulsingTemplate=1*0;

%LATER RUN USED FOR SECONDARY SCINTILLATION
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\18March_99_Ar_1_CF4\primary_230V\';
%Use fCF4 = 1.01 just to avoid ulterior condition in which the PM gain is modified due to small misadjustment of HVs.
%It actually replicates earlier values fairly well.
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%TakeAfterPulsingFromExt = 1; DIR_OUTPUTa = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\99_Ar_1_CF4\';

%%%%%%%%%%%%%%%%%
%%%    2.0%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\576V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\752V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%TakeAfterPulsingFromExt = 1; DIR_OUTPUTa = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\98_Ar_2_CF4\';    
%FIXME: this was apparently done automatically before. CHECK!

%%%%%%%%%%%%%%%%%
%%%    6.0%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\115V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\192V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\576V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\94_Ar_6_CF4\752V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;

%%%%%%%%%%%%%%%%%
%%%   10.0%  %%%% workspaces fully created
%%%%%%%%%%%%%%%%%

%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\0V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\7V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\19V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\38V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\76V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\114V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\191V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\306V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\383V\';
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\576V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
%DIR       = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\Argon_CF4\NewSetup\10bar\90_Ar_10_CF4\752V\';
%useAfterPulsingTemplate=1; storeAfterPulsingTemplate=1*0;
