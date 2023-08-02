%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% script for simulating OTPC data (DGD, 17/03/2019) %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% TO DO:                                                                  %
%                                                                         %
% -> Include a more realistic energy-loss                                 %
%                                                                         %
% TO IMPROVE:                                                             %
%                                                                         %
% -> Include Bethe-Bloch.                                                 %
% -> PMT time response totally approximate.                               %
% -> SinglePhoton width and noise from Q-analysis, not from A-analysis.   %
% -> Negative Tail of SinglePhoton probably unrealistic.                  %
% -> Include PMT-by-PMT variations.                                       % 
% -> Include global response parameters as workspaces and document, in    %
%    order to understand the various assumptions. All parameters          % 
%    should be revisited then.                                            % 
% -> Gap transit-time included as a convolution with a flat response.     %
%    Allow to change between exponential or flat, depending on the case   %
%    (CF4, noble)                                                         %
% -> Only direct light included.                                          %
% -> No CMOS offset. It requires maps...                                  %
% -> CMOS noise approximated through normal distribution                  %
% -> Include/introduce function for weighted histogram to include WePh    %
% electron by electron                                                    %
%                                                                         %
% 15 s with rebin 10, 10min with rebin 1 :(                               %
% Re-check CMOS numbers                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% I. INITIALIZATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OptFact = 10*10;                     %Correction factor to adjust optical gain

P      = 0.1;                     %#ok<GPFST> %[bar] Operational conditions under CF4 in early 2019 (P.A. tfg). 
Epart  = 5.485;                   % [MeV] Energy of alpha particles from 241Am (84.5% B.R. // 13% B.R., E = 5.442MeV).
Rpart  = 1.6/P * Epart / 5.4845;  % [cm]  Range of alpha particles for CF4, according to TRIM (Manuel Caamaño). FIXME.
Mpart  = 3727;                    % [MeV] 
z0     = 15;                      % [cm]  The maximum drift distance of the TPC is 15cm.

nsteps = round(Rpart/0.005);      % Sample initial ionization/scintillation over 50um steps.
theta0 = 45;                      % (deg). Default: 45
phi0   = 30;                      % (deg). Default: 30

eps    = 1e-4;                    % Just for numerics.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% II. GLOBAL RESPONSE PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% IONIZATION %%%%%
global Wi Ws tau1 tau2 Ratio12 FanoQ FanoS;
Wi      = 54;   % [eV] From "Properties of some gas mixtures used in tracking detectors", Archana Sharma, GSI Darmstadt, 1997.
Ws      = 500;  % [eV] At P~1bar, and E<60V/cm/bar. From Morozov's: https://doi.org/10.1016/j.nima.2010.07.001 , https://doi.org/10.1016/j.nimb.2010.01.012.
tau1    = 2;    % [ns] Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
tau2    = 40;   % [ns] Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
Ratio12 = 0.35; %      Morozov/Margato SECONDARY SCINTILLATION (only the dominant UV part -visible would be in between): doi: 10.1088/1748-0221/8/07/p07008.
FanoQ   = 0.2;  % Typical Fano factor for ionization    (unknown)
FanoS   = 0.2;  % Typical Fano factor for scintillation (unknown)

%%%%% DRIFT %%%%%
global P vd DL DT;
vd = 3.5;          % [cm/us]       Measured in Nausicaa1 at about 100 V/cm/bar (P. A. tfg), units changed below. 
DL = 250/sqrt(P);  % [um/sqrt(cm)] Simulated by Magboltz at about 100 V/cm/bar, units changed below.
DT = 250/sqrt(P);  % [um/sqrt(cm)] Simulated by Magboltz at about 100 V/cm/bar, units changed below.

%%%%% ANODE %%%%%
global hole pitch gap OptGain vdGap tau1Gap tau2Gap Ratio12Gap;

hole    = 0.05;       % [cm] (diameter), from S. Williams / R. de Oliveira. 
pitch   = 0.1;        % [cm] from S. Williams / R. de Oliveira (expected, but unconfirmed). 
gap     = 0.14;       % [cm] (thickness of the structure).
OptGain = 105*OptFact;% Experimentally obtained with meshes in early 2019 (P.A. tfg).
vdGap   = 10;         % [cm/us] preliminary simulated by Magboltz at about 10 kV/cm/bar (operating field in P.A. tfg was 25kV/cm/bar), units changed below.
tau1Gap = tau1;       % Assume same spectrum for primary and secondary scintillation.
tau2Gap = tau2;       % Assume same spectrum for primary and secondary scintillation.
Ratio12Gap = Ratio12; % Assume same spectrum for primary and secondary scintillation.

%%%%% PMT response %%%%%
global zPM rPM RPM QEPM Tmesh sigmaAovA1 sigmaNoise A1 sigmaTresponse meanTresponse Tbin BW;
zPM            = 24;        % [cm] PM z-position.
rPM            = 1.25;      % [cm] PM radius.
RPM            = 7.2;       % [cm] PM distance to center of flange.
QEPM           = 0.0885;    % Effective quantum efficiency for PMs (P.A. tfg).
Tmesh          = 0.71;      % Mesh optical transparency.
sigmaAovA1     = 1;         % From Q analysis (M.F. tfg). Approximate.  
sigmaNoise     = 0.2;       % In photon units (M.F. tfg). Approximate.
A1             = 5;         % [mV] Single photon amplitude. Approximate.
sigmaTresponse = 10;        % [ns] Width of PM time-response function. Approximate.
meanTresponse  = 20;        % [ns] Delay for PM time-response function. Approximate.
Tbin           = 10;        % [ns] Waveform sampling time.
BW             = 0.02;      % [GHz] Scope bandwidth for waveform digitization.
if(BW>0.1), BW = 0.1; end   % NOTE: it does not work for larger BW. This is not really a limitation, but it is strange.

%%%%% CMOS response %%%%%
global M lensN sigmaNph QECMOS rebin Npixel pixelSize;
M              = -1/23;     % magnification.
lensN          =  0.95;     % lens number.
sigmaNph       = 1.2;       % sigma in number of photons per pixel.
QECMOS         = 0.7*0.25;  % 70% of the 25% of visible light.
rebin          = 10;        % group pixels according to camera software.
Npixel         = 2048;      % number of pixel per line.
pixelSize      = 6.5e-4;    % pixel size [cm].

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% III. RUN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Adjust magnitudes of swarm parameters
DL    = DL*1e-4;      %cm/sqrt(cm)
DT    = DT*1e-4;      %cm/sqrt(cm)
vd    = vd*1e-3;      %cm/ns
vdGap = vdGap*1e-3;   %cm/ns

%Track generation
disp('track generation...')
[x, y, z, t, DeltaE]          = GenStraightTrack(z0, theta0, phi0, Rpart, Epart/Rpart, Mpart, nsteps);
disp('...done')
%Initial TPC response
disp('calculating primary TPC response...')
[xe, ye, ze, te, xPh, yPh, zPh, tPh, LPh] = TPCprimaryResponse(x, y, z, t, DeltaE);
disp('...done')
%TPC drift
disp('calculating TPC drift...')
[xeA, yeA, zeA, teA]          = TPCdrift(xe, ye, ze, te);
disp('...done')
%TPC anode
disp('calculating TPC anode response...')
[xePh, yePh, zePh, tePh, WePh] = TPCanode(xeA, yeA, zeA, teA);
disp('...done')
%PMT response
disp('calculating PMT response...')
[PML_wvf, PMR_wvf, PMU_wvf, PMD_wvf, T_wvf, tPML_S1, tPMR_S1, tPMU_S1, tPMD_S1,... 
    tPML_S2, tPMR_S2, tPMU_S2, tPMD_S2] = PMTresponse(xePh, yePh, zePh, tePh, WePh, xPh, yPh, zPh, tPh);
disp('...done')
disp('calculating CMOS response...')
NphXY = CMOSresponse(xePh, yePh, WePh);
disp('...done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% IV. DIAGNOSTIC PLOTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Draw ionization centers
figure; hold on;
subplot(2,2,1);
plot(x,   y, '.' ); hold on;
plot(xe, ye, 'or');
xaxis(0.9*min(x), 1.1*max(x) + eps); yaxis(0.9*min(y), 1.1*max(y) + eps); 
xlabel('x [cm]'); ylabel('y [cm]'); title('primary ionization');
legend('energy loss positions', 'e-ion pairs'); box; box;
subplot(2,2,2);
plot(x,   z,  '.'); hold on;
plot(xe, ze, 'or');
xaxis(0.9*min(x), 1.1*max(x) + eps); yaxis(0.9*min(z), 1.1*max(z) + eps);
xlabel('x [cm]'); ylabel('z [cm]'); title('primary ionization');
legend('energy loss positions', 'e-ion pairs'); box; box;
subplot(2,2,3);
plot(y,   z,  '.'); hold on;
plot(ye, ze, 'or');
xaxis(0.9*min(y), 1.1*max(y) + eps); yaxis(0.9*min(z), 1.1*max(z) + eps);
xlabel('y [cm]'); ylabel('z [cm]'); title('primary ionization');
legend('energy loss positions', 'e-ion pairs'); box; box;
subplot(2,2,4);
hist1D(te, 0:0.1:10);
xlabel('time[ns]'); ylabel('entries'); title('primary ionization');
box; box;

%Draw photon emission centers
figure; hold on;
subplot(2,2,1);
plot(x,   y, '.' ); hold on;
plot(xPh, yPh, 'sg');
xaxis(0.9*min(x), 1.1*max(x) + eps); yaxis(0.9*min(y), 1.1*max(y) + eps); 
xlabel('x [cm]'); ylabel('y [cm]'); title('primary scintillation');
legend('energy loss positions', 'photon emission points'); box; box;
subplot(2,2,2);
plot(x,   z,  '.'); hold on;
plot(xPh, zPh, 'sg');
xaxis(0.9*min(x), 1.1*max(x) + eps); yaxis(0.9*min(z), 1.1*max(z) + eps);
xlabel('x [cm]'); ylabel('z [cm]'); title('primary scintillation');
legend('energy loss positions', 'photon emission points'); box; box;
subplot(2,2,3);
plot(y,   z,  '.'); hold on;
plot(yPh, zPh, 'sg');
xaxis(0.9*min(y), 1.1*max(y) + eps); yaxis(0.9*min(z), 1.1*max(z) + eps);
xlabel('y [cm]'); ylabel('z [cm]'); title('primary scintillation');
legend('energy loss positions', 'photon emission points'); box; box;
subplot(2,2,4);
hist1D(tPh, 0:1:200); logy;
xlabel('time[ns]'); ylabel('entries'); title('primary scintillation');
box; box;

figure; hold on;
subplot(2,2,1);
plot(xeA, yeA, '.k');
xaxis(0.9*min(xeA), 1.1*max(xeA) + eps); yaxis(0.9*min(yeA), 1.1*max(yeA) + eps); 
xlabel('x_A [cm]'); ylabel('y_A [cm]');  title('primary ionization@anode');
box; box;
subplot(2,2,2);
plot(xeA, vd.*teA, '.k');
xaxis(0.9*min(xeA), 1.1*max(xeA) + eps); yaxis(0.9*min(ze), 1.1*max(ze) + eps);
xlabel('x_A [cm]'); ylabel('v_d \times t_A [cm]'); title('primary ionization@anode');
box; box;
subplot(2,2,3);
plot(yeA, vd.*teA, '.k');
xaxis(0.9*min(yeA), 1.1*max(yeA) + eps); yaxis(0.9*min(ze), 1.1*max(ze) + eps);
xlabel('y_A [cm]'); ylabel('v_d \times t_A [cm]'); title('primary ionization@anode');
box; box;
subplot(2,2,4);
hist1D(teA, 0:10:10000);
xlabel('time [ns]'); ylabel('entries'); title('primary ionization@anode');
box; box;

figure; hold on;
subplot(2,2,1);
plot(xePh, yePh, '.k');
xaxis(0.9*min(xePh), 1.1*max(xePh) + eps); yaxis(0.9*min(yePh), 1.1*max(yePh) + eps); 
xlabel('xePh [cm]'); ylabel('yePh [cm]');  title('secondary scintillation@anode');
box; box;
subplot(2,2,2);
plot(xePh, vd.*tePh, '.k');
xaxis(0.9*min(xePh), 1.1*max(xePh) + eps); yaxis(0.9*min(ze), 1.1*max(ze) + eps);
xlabel('xePh [cm]'); ylabel('v_d \times tePh [cm]'); title('secondary scintillation@anode');
box; box;
subplot(2,2,3);
plot(yePh, vd.*tePh, '.k');
xaxis(0.9*min(yePh), 1.1*max(yePh) + eps); yaxis(0.9*min(ze), 1.1*max(ze) + eps);
xlabel('yePh [cm]'); ylabel('v_d \times tePh [cm]'); title('secondary scintillation@anode');
box; box;
subplot(2,2,4);
hist1D(tePh, 0:10:10000);
xlabel('tePh [ns]'); ylabel('entries'); title('secondary scintillation@anode');
box; box;

figure; hold on;
subplot(2,2,1)
plot(T_wvf, PML_wvf, '-'); xaxis(-1000, max(T_wvf)); yaxis(min(PML_wvf), max(PML_wvf)); xlabel('time [ns]'); ylabel('V [mV]'); title('PM-left');
subplot(2,2,2)
plot(T_wvf, PMR_wvf, '-'); xaxis(-1000, max(T_wvf)); yaxis(min(PMR_wvf), max(PMR_wvf)); xlabel('time [ns]'); ylabel('V [mV]'); title('PM-right');
subplot(2,2,3)
plot(T_wvf, PMU_wvf, '-'); xaxis(-1000, max(T_wvf)); yaxis(min(PMU_wvf), max(PMU_wvf)); xlabel('time [ns]'); ylabel('V [mV]'); title('PM-up');
subplot(2,2,4)
plot(T_wvf, PMD_wvf, '-'); xaxis(-1000, max(T_wvf)); yaxis(min(PMD_wvf), max(PMD_wvf)); xlabel('time [ns]'); ylabel('V [mV]'); title('PM-down');

sizeNphXY  = size(NphXY); 
NpixXrebin = 1:sizeNphXY(1);
NpixYrebin = 1:sizeNphXY(2); 
[NpixX_mesh, NpixY_mesh] = meshgrid(NpixXrebin, NpixYrebin);
figure;
surf_D(NpixX_mesh, NpixY_mesh, NphXY); xlabel('pixel x'); ylabel('pixel y'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% V. SUMMARY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp(['<n_e>             = ',  num2str(round(Epart*1e+6/Wi))]);
disp(['n_e               = ',  num2str(length(xe))]);
disp(['<n_{ph,S1}>       = ',  num2str(round(Epart*1e+6/Ws))]);
disp(['n_{ph,S1}         = ',  num2str(length(xPh))]);

disp(['<n_{ph,S2}>       = ',  num2str(round(sum(WePh)))]);
disp(['m_\gamma          = ',  num2str(sum(WePh)/length(xe))]);

disp(['<n_{det,S1}/PM>   = ',  num2str(sum([length(tPML_S1),length(tPMR_S1),length(tPMU_S1),length(tPMD_S1)])/4)]);
disp(['<n_{det,S2}/PM>   = ',  num2str(sum([length(tPML_S2),length(tPMR_S2),length(tPMU_S2),length(tPMD_S2)])/4)]);
disp(['<n_{det,S2}/e/PM> = ',  num2str(sum([length(tPML_S2),length(tPMR_S2),length(tPMU_S2),length(tPMD_S2)])/4/length(xe))]);

toc; mosaic;

return;
