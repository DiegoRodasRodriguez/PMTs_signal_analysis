%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% initialization file for ana_FATGEM_DEV %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Ndata Qdata sigmaNdata;         %#ok<NUSED> %For fit

%Qsp = 1/3.3; %To go from 1.35kV to 1.15kV (unused explicitely in code)

Qth = [];
Ath = [];

%%%%% 2021 run %%%%%
%FAT-GEM, EL-SCAN
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\2500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 100;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\3000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 200;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\3500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 300;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\4000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 400;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\4500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
end

%FAT-GEM, DRIFT-SCAN
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\200_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4.2; TevtMIN = 3.7; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\300_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4.1; TevtMIN = 3.7; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\520_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4; TevtMIN = 3.5; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\620_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4; TevtMIN = 3.5; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\620_V_redo\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4; TevtMIN = 3.5; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\800_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 4; TevtMIN = 3.5; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\1000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1; Qth = 500;
end

%PEN-GEM(A), EL-SCAN
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\2000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 200;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\2500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.07;  TfevtMIN = 0.07;  Qth = 70;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\3000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 100;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\3500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 200;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\4000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 200;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\4400_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 1.5; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 200;
end

if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Argon\4bar_13May\EL_3000V\4000_V\'))
     stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 100;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Argon\4bar_13May\EL_6000V\7500_V\'))
    stdMIN     = 0.28; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = inf; AMAX = inf; TevtMAX = 2.2; TevtMIN = 1.6; TrevtMIN = 0.2;  TfevtMIN = 0.2;  Qth = 0;
end

%PEN-GEM(B), EL-SCAN
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\2000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 4; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.05;  TfevtMIN = 0.05;  Qth = 100;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\2500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 4; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.05; TfevtMIN = 0.05; Qth = 70;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\3000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 4; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 100;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\3500_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 4; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 150;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\4000_V\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 4; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.1;  TfevtMIN = 0.1;  Qth = 150;
end

%CHRISTMAS 2021, FIRST STRUCTURE (canonical, check name)
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe_ChristmasCampaign\MCA_to_compare\Xenon\150Vcmbar5.0kVcmbar\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.06;  TfevtMIN = 0.06;
end

%FAT-GEM++ (for reference)
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF65_EL4.35\Filter\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.06;  TfevtMIN = 0.06;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF65_EL4.35\NO_Filter\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.06;  TfevtMIN = 0.06;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF167_EL4.35\NO_Filter\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.06;  TfevtMIN = 0.06;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF167_EL4.35\Filter\'))
    stdMIN     = 0; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = inf; TevtMIN = 0; TrevtMIN = 0.06;  TfevtMIN = 0.06;
end
%FAT-GEM++ (6bar)
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_25\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_28\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_315\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_34\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_37\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
end
%FAT-GEM++ (8bar)
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_173\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_19\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_21\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_23\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_25\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_28\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_315\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_34\'))
    stdMIN     = 0.5; stdMAX     = inf; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 4.5; TevtMIN = 3.2; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0; 
end
%FAT-GEM++ (4bar)
if(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\PEDOT_REFL\2_5_EL_50_DR_1150_PMT\'))
    stdMIN     = 0.5; stdMAX     = 1.5; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 11; TevtMIN = 9; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0;
elseif(strcmp(DIR, 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\PEDOT_REFL\3_15_EL_50_DR_1150_PMT\'))
   stdMIN     = 0.5; stdMAX     = 1.5; BasevtMAX     = inf; BasevtMAX_HC  = inf; AssyevtMAX    = inf; AssyevtMAX_HC = 5; AMAX = inf; TevtMAX = 11; TevtMIN = 9; TrevtMIN = 0.2;  TfevtMIN = 0.2; Qth = 100*0; Ath = 4*0;
end
