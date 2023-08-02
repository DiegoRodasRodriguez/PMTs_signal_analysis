%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Initialization script for analyzing 1PMT FATGEM-data from       %%%%%
%%%%% LED, with automatic fit (companion of AnaFATGEM)                %%%%%
%%%%%                                                                 %%%%%
%%%%%                   (DGD 31/Aug/2020)                             %%%%%
%%%%%               (final version 10/Dec/2022)                       %%%%%
%%%%%                                                                 %%%%%
%%%%%   USE NOTES:                                                    %%%%%
%%%%%                                                                 %%%%%
%%%%%   - Initialization parameters saved in this file for later      %%%%%
%%%%%   reanalysis.                                                   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%INCLUDES FAT & PEN-GEM RUNS FROM 2021

% 2021
%%%%%%%%%%%%%%%%%
% FAT-GEM with anode-mesh, as in 2019 but with larger metalization; P = 2bar in xenon (smart trigger on signal window not yet used)
%%%%%%%%%%%%%%%%%

%isTriggerSignalIn1 = 0;
%nCh        = 1;

%EL-SCAN
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\2500_V\';  %%%Eres=38.5% Q=290
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\3000_V\';   %%%Eres=37% Q=490
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\3500_V\';  %%%Eres=32% Q=690
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\4000_V\';  %%%Eres=27% Q=930
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\DF_620V\4500_V\';  %%%Eres=27% Q=1170

%DRIFT-SCAN
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\200_V\';       %%%Eres=27.5% Q=1630
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\300_V\';       %%%Eres=23.4% Q=1770
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\520_V\';       %%%Eres=22.7% Q=1690
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\620_V\';       %%%Eres=25.4% Q=1630
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\620_V_redo\';  %%%Eres=25.4% Q=1330
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\800_V\';       %%%Eres=26.1% Q=1310
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\FAT_GEMS\Xenon\2bar_10May\EL_5050V\1000_V\';       %%%Eres=31.3% Q=1190

%%%%%%%%%%%%%%%%%
% PEN-GEM(A) with double mesh; P = 2bar in xenon (smart trigger on signal window, trigger signal stored in first channel)
%%%%%%%%%%%%%%%%%
%isTriggerSignalIn1 = 1;

%nCh          = 2;
%EL-SCAN
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\2000_V\'; %%% TOO LOW
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\2500_V\';    %%%Eres = 0.69, Q=129
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\3000_V\';   %%%Eres = 0.48, Q=214
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\3500_V\';   %%%Eres = 0.45, Q=305
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\4000_V\';   %%%Eres = 0.43, Q=384
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\DF_1000V_redo\4400_V\';  %%%Eres = 0.47, Q=371

%DRIFT-SCAN
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\200_V\';
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\400_V\';
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\600_V\';
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\800_V\';
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\1000_V\';
%DIR           = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_13May\EL_4000V\1000_V_final\';

% PEN-GEM(A) with double mesh; P = 4bar in argon (smart trigger on signal window, trigger signal stored in first channel)

%EL 3000V, Cathode 4000V
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Argon\4bar_13May\EL_3000V\4000_V\';
%EL 6000V, Cathode 7500V
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Argon\4bar_13May\EL_6000V\7500_V\';

%%%%%%%%%%%%%%%%%%%
% PEN-GEM(B) with anode-mesh; P = 2bar in xenon (smart trigger on signal window, trigger signal stored in first channel)
%%%%%%%%%%%%%%%%%%%
%isTriggerSignalIn1 = 1;
%nCh        = 2;
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\2000_V\';  %%% TOO LOW
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\2500_V\';  %%%Eres = 0.43, Q=111
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\3000_V\';  %%%Eres = 0.38, Q=160.7
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\3500_V\';  %%%Eres = 0.36, Q=218.4
%DIR            = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe\PEN_GEMS\Xenon\2bar_14May_coating\DF_1000V\4000_V\';  %%%Eres = 0.32, Q=274

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Christmas 2021 campaign %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% isTriggerSignalIn1 = 1;
% nCh        = 2;
% DIR        = 'E:\HOME_RareEventsGroup\SETUPS\FALCON\DATA\2021Xe_ChristmasCampaign\MCA_to_compare\Xenon\150Vcmbar5.0kVcmbar\';
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% FAT-GEM++ 2022 campaign %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%P=2
% Data taken for comparison with MCA and relative normalization (looks like no data: Check again)
%isTriggerSignalIn1 = 1;
%nCh        = 2;
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF65_EL4.35\Filter\';
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF65_EL4.35\NO_Filter\';
%DIR        = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF167_EL4.35\NO_Filter\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\Card_to_compare_with_MCA\Xenon\2bar\DF167_EL4.35\Filter\';

%P=6
%isTriggerSignalIn1 = 0;
%nCh         = 1;
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_2\'; Only 5kevts saved, by mistake
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_25\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_28\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_315\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_34\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\6bar_Dec22_2021\DF_100\EL_37\';

%P=8
%isTriggerSignalIn1 = 0;
%nCh         = 1;
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_173\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_19\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_21\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_23\';
DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_25\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_28\';   %bad fit
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_315\';
%DIR         = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\calibration\8bar_Dec22_2022\DF_75\EL_34\';

%P=4 %(Clevios run)
%isTriggerSignalIn1 = 0;
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\PEDOT_REFL\2_5_EL_50_DR_1150_PMT\';
%DIR          = 'E:\HOME_RareEventsGroup\SETUPS\SWAN\DATA\Falcon_fatgems\PEDOT_REFL\3_15_EL_50_DR_1150_PMT\';