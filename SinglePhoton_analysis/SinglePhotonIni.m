%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialization script for analyzing single-photon PMT data from %%%%%
%%%%% LED, with automatic fit (companion of AnaSinglePhoton)          %%%%%
%%%%%                                                                 %%%%%
%%%%%                   (DGD 31/Aug/2020)                             %%%%%
%%%%%               (final version 10/Dec/2022)                       %%%%%
%%%%%                                                                 %%%%%
%%%%%   USE NOTES:                                                    %%%%%
%%%%%                                                                 %%%%%
%%%%%   - Initialization parameters saved in this file for later      %%%%%
%%%%%   reanalysis.                                                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%--------------------------PMT-Falcon (2018)---------------------------
%----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Reference data from Manuel Fontaíña with PMT7378,        %%%%%%%
%%%%%%%%%%%%%%%%%%%% NEXT bases and HV = 1350V %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DataType = 'OSCI';
%nCh = 1; isLEDsignal1 = 0;   %no LED signal stored (but triggered on it)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\SinglePhoton\LEDscan&SPsearch\';
%FILE = '075mv60kevt';  Nfiles = 6; %Def: 6. Available up to 6
% FILE = '15mv100kevt';  Nfiles = 1; %Def: 7. Available up to 9 (last ones seem to drift)
% FILE = '3mv60kevt';    Nfiles = 6;  %Def: 6. Available up to 6
%AmplitudeWindow = [481,500];         %Window to look for signal amplitude
%ChargeWindow    = [470,540];         %Window to look for signal charge


%----------------------------------------------------------------------
%---------------------------PMT1-SWAN (2020)---------------------------
%----------------------------------------------------------------------

%%%%%% PMT1:RUN1 (PMT7378@1350 DS-bases, swanXePaper2020)  %%%%%%%%%%%%
%decent: old pulser, some x-talk
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;   %default, LED is second channel
%PMT1
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_1_YELLOW\75\';      %Does not converge
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_1_YELLOW\50\';      %1.I %P0 = 41%     // <Q> = 0.886pC  // sigmaQ = 0.565pC %g// this looks quite good
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_1_YELLOW\20\';      %1.II P0 = 15.6%   // <Q> = 1.158pC  // sigmaQ = 0.38pC  %g//
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_1_YELLOW\5\';       %1.III P0 = 15.8%  // <Q> = 1.163pC  // sigmaQ = 0.46pC  %g//
%FILE = 'wave';
%AmplitudeWindow = [670,740];           %Window to look for signal amplitude
%ChargeWindow    = [670,740];           %Window to look for signal charge
%
%<Q>_1 at 1.35kV  ==  1.02 +- 0.19  (take only first two, sigma == rms)
%G = G0*exp(7.214*(V-V0))
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
%%%%%% PMT1:RUN2 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%bad: old pulser, too high cross-talk
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\YellowPMT1\30percent\';  %1.IV bad %P0 = 44%  // <Q> = 1.46pC   // sigmaQ = 0.10pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\YellowPMT1\50percent\';  %1.V bad %P0 = 76%   // <Q> = 1.629pC  // sigmaQ = 0.10pC
%FILE = 'wave';
%AmplitudeWindow = [3930,3980];           %Window to look for signal amplitude
%ChargeWindow    = [3930,3980];           %Window to look for signal charge
%
%%%%%% PMT1:RUN3 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%bad: new pulser, too high cross-talk (in disagreement with previous value around 1pC)
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;   %default, LED is first channel
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow\30percent\';  %1.VI   P0 = 64.5%  // <Q> = 0.547pC   // sigmaQ = 0.346pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow\50percent\';  %1.VII  P0 = 62%    // <Q> = 0.497pC   // sigmaQ = 0.320pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow\70percent\';  %1.VIII P0 = 77%    // <Q> = 0.57pC    // sigmaQ = 0.31pC
%FILE = 'wave';
%AmplitudeWindow = [380,500];           %Window to look for signal amplitude
%ChargeWindow    = [380,500];           %Window to look for signal charge
%
%%%%%% PMT1:RUN4 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow_2\30percent\';  %1.IX P0 = 50%    // <Q> = 0.876pC   // sigmaQ = 0.659pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow_2\50percent\';  %1.X  P0 = 67%    // <Q> = 0.696pC   // sigmaQ = 0.604pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_1_Yellow_2\70percent\';  %1.XI P0 = 84%    // <Q> = 0.68pC    // sigmaQ = 0.58pC
%FILE = 'wave';
%AmplitudeWindow = [380,550];           %Window to look for signal amplitude
%ChargeWindow    = [380,550];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT 1 %%%%%%%%%%%%%%%%%%%%%%%%
%
% Previous:
% <Q>_1 at 1.35kV  ==  1.02 +- 0.19  (take only first two, sigma == rms)
%
% New (take best fits 1.I and 1.IX and assign error to be compatible with previous result at 2sigma:
% <Q>_1 at 1.35kV  ==  0.88 +- 0.07
% G = G0*exp(7.214*(V-V0))
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: Spectrum 1.IX
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%---------------------------PMT2-SWAN (2020)---------------------------
%----------------------------------------------------------------------

%%%%%% PMT2:RUN1 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;
%PMT2
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_2_BLUE\75\'; %2.I   P0 = 82%    // <Q> = 0.467pC   // sigmaQ = 0.373pC %g// this looks beautiful!
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_2_BLUE\50\'; %2.II  P0 = 35.7%  // <Q> = 0.548pC   // sigmaQ = 0.428pC %g//
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_2_BLUE\20\'; %2.III P0 = 10%    // <Q> = 0.758pC   // sigmaQ = 0.4pC %g-b// large number of photons
%FILE = 'wave';
%AmplitudeWindow = [670,740];           %Window to look for signal amplitude
%ChargeWindow    = [670,740];           %Window to look for signal charge
%
%<Q>_2 at 1.35kV  ==  0.51 +- 0.06  (take only first two, sigma == rms)
%G = G0*exp(7.124*(V-V0))
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
%%%%%% PMT2:RUN2 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2; %DOUBLE-CHECK!!
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\BluePMT2\30percent\';  %2.IV P0 = 6%   // <Q> = 0.579pC   // sigmaQ = 0.10908 // too many photons-bad fit (too narrow Q-dist)
%DIR ='E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\BluePMT2\50percent\';   %2.V  P0 = 6%   // <Q> = 0.526pC   // sigmaQ = 0.112pC // too many photons-bad fit (too narrow Q-dist)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\BluePMT2\70percent\';  %2.VI P0 = 80%  // <Q> = 0.62pC    // sigmaQ = 0.32pC  //
%FILE = 'wave';
%AmplitudeWindow = [3930,3970];           %Window to look for signal amplitude
%ChargeWindow    = [3930,3970];           %Window to look for signal charge
%
%%%%%% PMT2:RUN3 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_2_Blue\30percent\';  %2.VII  P0 = 57%    // <Q> = 0.448pC   // sigmaQ = 0.151pC -->Peak changes with time above 10kevts!
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_2_Blue\50percent\';  %2.VIII P0 = 69%    // <Q> = 0.422pC   // sigmaQ = 0.27pC good!
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_2_Blue\70percent\';  %2.IX   P0 = 82%    // <Q> = 0.506pC   // sigmaQ = 0.1pC too narrow!
%FILE = 'wave';
%AmplitudeWindow = [380,500];           %Window to look for signal amplitude
%ChargeWindow    = [380,500];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT 2 %%%%%%%%%%%%%%%%%%%%%%%%
%
% Previous:
%<Q>_2 at 1.35kV  ==  0.51 +- 0.06  (take only first two, sigma == rms)
%
% New (take best fits 2.I and 2.VIII and give rms to include previous value at 2sigma)
%<Q>_2 at 1.35kV  ==  0.45 +- 0.03
%G = G0*exp(-7.124*0.2)
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: Spectrum 2.VIII
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%---------------------------PMT3-SWAN (2020)---------------------------
%----------------------------------------------------------------------

%%%%%% PMT3:RUN1 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_3_RED\75\'; %3.I   P0 = 89.9%    // <Q> = 0.858pC   // sigmaQ = 0.528pC %g// this looks beautiful!
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_3_RED\50\'; %3.II  P0 = 61.2%    // <Q> = 0.940pC   // sigmaQ = 0.761pC %g// quite nice too
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_3_RED\20\'; %3.III P0 = 31.8%    // <Q> = 1.019pC   // sigmaQ = 0.819pC %g// quite nice too
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_3_RED\5\';  %3.IV  P0 = 12.2%    // <Q> = 1.07pC    // sigmaQ = 0.73pC  %g// decent
%FILE = 'wave';
%AmplitudeWindow = [670,740];           %Window to look for signal amplitude
%ChargeWindow    = [670,740];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT 3 %%%%%%%%%%%%%%%%%%%%%%%%
%
% Previous == New (no new fit, as PM got replaced):
%<Q>_3 at 1.35kV  ==  0.97 +- 0.09  (take the four, sigma == rms)
%G = G0*exp(6.966*(V-V0))
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: Spectrum 3.II
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%--------------------------PMT3b(5070)-SWAN (2020)---------------------
%----------------------------------------------------------------------

%%%%%% PMT3b:RUN1 (PMT5070@1350 DS-bases, swanArCF4Paper2020)CHECK!! %%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\RedPMT3\30percent\';  %3b.I    P0 = 6.8% // <Q> = 0.576pC  // sigmaQ = 0.583pC (too many photons)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\RedPMT3\50percent\';  %3b.II   P0 = 13%  // <Q> = 0.589pC  // sigmaQ = 0.335pC (too many photons -unstable against integration window)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Apr2021\NewSetup\RedPMT3\70percent\';  %3b.III  P0 = 33%  // <Q> = 0.683pC  // sigmaQ = 0.453pC (too many photons -unstable against integration window)
%FILE = 'wave';
%AmplitudeWindow = [3930,3970];           %Window to look for signal amplitude
%ChargeWindow    = [3930,3970];           %Window to look for signal charge
%
%%%%%% PMT3b:RUN2 (PMT5070@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_3_Red\20percent\';   %3b.IV  P0 = 61%    // <Q> = 0.443pC   // sigmaQ = 0.376pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_3_Red\50percent\';   %3b.V   P0 = 68%    // <Q> = 0.417pC   // sigmaQ = 0.345pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_3_Red\50percent_2\'; %3b.VI  P0 = 70%    // <Q> = 0.498pC   // sigmaQ = 0.344pC
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_3_Red\80percent\';   %3b.VII P0 = 93%    // <Q> = 0.327pC   // sigmaQ = 0.273pC //(too little photons
%FILE = 'wave';
%AmplitudeWindow = [380,500];           %Window to look for signal amplitude
%ChargeWindow    = [380,500];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT 3b %%%%%%%%%%%%%%%%%%%%%%%%
%
% New (no previous fit, as PM got replaced):
%<Q>_3 at 1.35kV  ==  0.45 +- 0.04  (take the first 3, sigma == rms)
%G = G0*exp(-7.1*0.2);      %Average trend from the others
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: Spectrum 3b.IV
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%--------------------------PMT4-SWAN (2020)----------------------------
%----------------------------------------------------------------------

%%%%%% PMT4:RUN1 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 2;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\3_take\75\'; %too many zeros
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\3_take\50\'; %too many zeros
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\3_take\30\'; %4.I   P0 = 4%     // <Q> = 0.614pC   // sigmaQ = 0.215pC        %b// little zeros
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\2_take\75\'; %4.II  P0 = 20%    // <Q> = 0.406pC   // sigmaQ = 0.231pC        %g// zeros reasonably well catched
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\2_take\50\'; %4.III P0 = 14%    // <Q> = 0.435pC   // sigmaQ = 0.1pC (forced) %b// too narrow response
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\2_take\30\'; %4.IV  P0 =  5%    // <Q> = 0.608pC   // sigmaQ = 0.1pC (forced) %b// too narrow response
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\calibration\PMT_Volt_1350V\PMT_4_GREEN\1_take\50\'; %4.V   P0 = 18%    // <Q> = 0.486pC   // sigmaQ = 0.27pC         %b// pedestal a bit high up
%FILE = 'wave';
%AmplitudeWindow = [670,740];           %Window to look for signal amplitude
%ChargeWindow    = [670,740];           %Window to look for signal charge
%
%<Q>_4 at 1.35kV  ==  0.50 +- 0.09  (take all, sigma == rms)
%G = G0*exp(-7.1*0.2); [discard strange measurement yielding G0*exp(5.093*(V-V0))]
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
%%%%%% PMT4:RUN1 (PMT7378@1350 DS-bases, swanXePaper2020) %%%%%%%%%%%%%
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_4_Green\30percent\';  %4.VI   P0 = 65%    // <Q> = 0.44pC   // sigmaQ = 0.26pC  // best so far
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_4_Green\50percent\';  %4.VII  P0 = 85%    // <Q> = 0.57pC   // sigmaQ = 0.11pC  // too narrow response
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\2020run\PMT_Calibration_Jun2021\NewSetup\PMT_4_Green\70percent\';  %4.VIII P0 = 89%    // <Q> = 0.40pC   // sigmaQ = 0.267pC //
%AmplitudeWindow = [380,500];           %Window to look for signal amplitude
%ChargeWindow    = [380,500];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT 4 %%%%%%%%%%%%%%%%%%%%%%%%
%
% Previous:
%<Q>_4 at 1.35kV  ==  0.50 +- 0.09
%
% New (take best fits 4.II, 4.VI and 4.VIII and give rms to include previous value at 2sigma)
%<Q>_4 at 1.35kV  ==  0.42 +- 0.04
%G = G0*exp(-7.1*0.2);
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: Spectrum 4.VI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------
%--------------------------PMT-Falcon (2021)---------------------------
%----------------------------------------------------------------------

%%%%%% Different from 2018 , quite probably the same as 2022, RE-CHECK %%
%%%PMT-FALCON (May 2021)  -> SP peak easily visible in amplitude!. Consider fitting it for double-checking zeros. Too high x-talk blurries charge estimate though.
%DataType = 'CAEN';
%nCh = 2; isLEDsignal1 = 1;
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Calibration_pmt_falcon\20percent\'; %F.I   P0 = 27%    // <Q> = 0.99pC   // sigmaQ = 0.66pC  %g// reasonable (x-talk a bit too high)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Calibration_pmt_falcon\50percent\'; %F.II  P0 = 76%    // <Q> = 1.228pC  // sigmaQ = 0.370pC %g// looks good (x-talk a bit too high)
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Calibration_pmt_falcon\80percent\'; %F.III P0 = 70%    // <Q> = 0.959pC  // sigmaQ = 0.605pC %g// looks good (x-talk a bit too high)
%FILE = 'wave';
%AmplitudeWindow = [1950,2000];           %Window to look for signal amplitude
%ChargeWindow    = [1950,2000];           %Window to look for signal charge
%
%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT-F %%%%%%%%%%%%%%%%%%%%%%%%
%
% Previous:
%<Q>_F at 1.35kV  ==  0.60 +- 0.21
%
% New (take mean and std to be compatible with previous value at 2sigma)
%<Q>_F at 1.35kV  ==  1.06 +- 0.20
%G = G0*exp(-7.1*0.2);
%(note: gain is better parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty. This would give a better agreement with Xe paper by 20%)
%
% FOR PAPER: calibrate again trying to avoid x-talk problem and use
% amplitude analysis to constrain parameters of charge fit.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------
%--------------------------PMT-Falcon (November-2022)------------------
%----------------------------------------------------------------------

%%%%%% Improved over 2022 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DataType = 'CAEN';
nCh = 2; isLEDsignal1 = 1;

DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration_Falcon_05Dec2022\80sig_20zero\'; % P0 = 67%    // <Q> = 1.03pC   // sigma_Q/<Q> = 0.71  // very good
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration_Falcon_05Dec2022\50sig_50zero\';  % P0 = 79%    // <Q> = 0.94pC   // sigma_Q/<Q> = 0.77  // very good
%DIR = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration_Falcon_05Dec2022\20sig_80zero\'; % P0 = 90%    // <Q> = 0.70pC   // sigma_Q/<Q> = 1.02  // good (but not considered in fit)

FILE = 'wave';
AmplitudeWindow = [700,780];           %Window to look for signal amplitude
ChargeWindow    = [700,780];           %Window to look for signal charge

%%%%%%%%%%%%%%%%%%%%%%%%%% FINAL RESULTS PMT-F %%%%%%%%%%%%%%%%%%%%%%%%
%
% New
%<Q>_F at 1.35kV  ==  0.985 +- 0.06
%(gain is well parameterized by a power-law, yielding a factor 3.3 over 200V, with a 10% uncertainty)
% Gain = k * V^7.63
%
% FOR PAPER: calibrate again trying to avoid x-talk problem and use
% amplitude analysis to constrain parameters of charge fit.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%trigTime        = (AmplitudeWindow(1) + AmplitudeWindow(2))/2 - 100;
