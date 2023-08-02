clear all
close all

disp('                                                             ')
disp('           *------------------------------------------------*')
disp('           *    Macro for diagnose signals from scope       *')
disp('           *------------------------------------------------*')
disp('                                                             ')
pause(1)

%HERE GIVE FILE NAME

%-------------------------analysis for RPC FEE
base_name     = 'F:\DANI_ELE_ANALISIS\';
%file_name     = 'CosmicRays_Avalanch_LogicTrigg_081113';
file_name     = 'CosmicRays_LogicTrigg_081106';
extension     = '.mat';
file          = [base_name file_name extension];
eval(['load ', file]);

%
N_points     = 500; %Points per waveform

N_first      =   1;
N_plots      =  10;
N_waves_plot = 100;

ndown = (N_first-1)*N_points + 1;

for i=1:N_plots

    ndown = (N_first-1)*N_points + 1 + (i-1)*N_points*N_waves_plot;
    nup   = ndown + N_points*N_waves_plot - 1;
    
    figure;
    plot(ndown:1:nup,U(ndown:nup,3),'-'); hold on;
    plot(ndown:1:nup,U(ndown:nup,4),'g-');
    plot(ndown:1:nup,U(ndown:nup,5),'r-');
    legend('LVDS','RPC','scintillator');
end

mosaic

return;