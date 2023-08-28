%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% script to simulate PMT reponse (DGD 18/03/19) %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  function [PML_wvf, PMR_wvf, PMU_wvf, PMD_wvf, Twvf, tPML_S1, tPMR_S1,...
%    tPMU_S1, tPMD_S1, tPML_S2, tPMR_S2, tPMU_S2, tPMD_S2] = ...          %
%    PMTresponse(xS2, yS2, zS2, tS2, WS2, xS1, yS1, zS1, tS1)             %   
%                                                                         %
% INPUT:                                                                  %
% xS2, yS2, zS2, tS2, WS2: position time and weight for S2                %
% xS1, yS1, zS1, tS1:      position time and weight for S1                %
%                                                                         %
% OUTPUT:                                                                 %
% PML,R,U,D_wvf: PM wvf's                                                 %
% tPML,R,U,D_S1: S1 PM times PMR_wvf, PMU_wvf, PMD_wvf,                   %
% tPML,R,U,D_S2: S2 PM times PMR_wvf, PMU_wvf, PMD_wvf,                   %
%                                                                         %
%                                                                         %
% PARAMETERS:                                                             %
% zPM            : z distance to PM plane    [cm]                         %
% rPM            : radius of PM photocathode [cm]                         %
% RPM            : radial distance to PM center in PM plane [cm]          %
% QEPM           : effective quantum efficiency                           %
% Tmesh          : mesh transparency                                      %      
% sigmaAovA1     : single-photon relative width                           %
% sigmaNoise     : sigma of noise in photon units                         %
% A1             : conversion from photons to amplitude                   %
% sigmaTresponse : single-photon time response function width [ns]        %
% Tresponse      : single-photon time response function delay [ns]        %
% Tbin           : time bin [ns]                                          %
% BW             : acquisition bandwidth [GHz]                            %
%                                                                         %
% gap            : size of S2 gap [cm]                                    %
% vdGap          : electron velocity in S2 gap [cm/ns]                    %
%                                                                         %
% NOTES/TO DOs:                                                           %
% * Gap transit time introduced through a convolution with a flat         %
% distribution. True for EL. For CF4 an exponential will be more suited   %
% Include this as a parameter.                                            %
% * No reflections included, only direct light.                           %
%                                                                         %  
% REFERENCE SYSTEM                                                        %
%                                                   ^  y                  %
%              o                (UP)                |                     %
%     LEFT)  o   o   (RIGHT)                x   <---|                     %
%              o                (DOWN)                                    %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [PML_wvf, PMR_wvf, PMU_wvf, PMD_wvf, Twvf, tPML_S1, tPMR_S1,...
    tPMU_S1, tPMD_S1, tPML_S2, tPMR_S2, tPMU_S2, tPMD_S2] = ...
    PMTresponse(xS2, yS2, zS2, tS2, WS2, xS1, yS1, zS1, tS1)

global zPM rPM RPM QEPM Tmesh sigmaAovA1 sigmaNoise A1 sigmaTresponse meanTresponse Tbin BW;
global gap vdGap;

clight  = 30;                  % [cm/ns]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% GENERATE TIME ENTRIES FOR DETECTED PHOTONS %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Primary scintillation weights/probabilities per time entry
WS1L      = zeros(size(tS1)); WS1R      = zeros(size(tS1)); WS1U      = zeros(size(tS1)); WS1D      = zeros(size(tS1));
lengthS1L = zeros(size(tS1)); lengthS1R = zeros(size(tS1)); lengthS1U = zeros(size(tS1)); lengthS1D = zeros(size(tS1)); 

for i=1:length(tS1)
    lengthS1L(i) = sqrt((zS1(i)-zPM)^2 + (xS1(i) - RPM)^2 +  yS1(i)^2);
    lengthS1R(i) = sqrt((zS1(i)-zPM)^2 + (xS1(i) + RPM)^2 +  yS1(i)^2);
    lengthS1U(i) = sqrt((zS1(i)-zPM)^2 +  xS1(i)^2        + (yS1(i) - RPM)^2); 
    lengthS1D(i) = sqrt((zS1(i)-zPM)^2 +  xS1(i)^2        + (yS1(i) + RPM)^2);    
    WS1L(i)      = 1/4 * rPM^2/lengthS1L(i)^2;
    WS1R(i)      = 1/4 * rPM^2/lengthS1R(i)^2;    
    WS1U(i)      = 1/4 * rPM^2/lengthS1U(i)^2;    
    WS1D(i)      = 1/4 * rPM^2/lengthS1D(i)^2;
end    
WS1L = WS1L * QEPM * Tmesh;
WS1R = WS1R * QEPM * Tmesh;
WS1U = WS1U * QEPM * Tmesh;
WS1D = WS1D * QEPM * Tmesh;

%Secondary scintillation weights/probabilities per time entry
WS2L      = zeros(size(tS2)); WS2R      = zeros(size(tS2)); WS2U      = zeros(size(tS2)); WS2D      = zeros(size(tS2));
lengthS2L = zeros(size(tS2)); lengthS2R = zeros(size(tS2)); lengthS2U = zeros(size(tS2)); lengthS2D = zeros(size(tS2)); 

for i=1:length(tS2)
    lengthS2L(i) = sqrt((zS2(i)-zPM)^2 + (xS2(i) - RPM)^2 +  yS2(i)^2);
    lengthS2R(i) = sqrt((zS2(i)-zPM)^2 + (xS2(i) + RPM)^2 +  yS2(i)^2);
    lengthS2U(i) = sqrt((zS2(i)-zPM)^2 +  xS2(i)^2        + (yS2(i) - RPM)^2); 
    lengthS2D(i) = sqrt((zS2(i)-zPM)^2 +  xS2(i)^2        + (yS2(i) + RPM)^2);    
    WS2L(i)      = WS2(i) * 1/4 * rPM^2/lengthS2L(i)^2;
    WS2R(i)      = WS2(i) * 1/4 * rPM^2/lengthS2R(i)^2;    
    WS2U(i)      = WS2(i) * 1/4 * rPM^2/lengthS2U(i)^2;    
    WS2D(i)      = WS2(i) * 1/4 * rPM^2/lengthS2D(i)^2;
end    
WS2L = WS2L * QEPM * Tmesh^2;
WS2R = WS2R * QEPM * Tmesh^2;
WS2U = WS2U * QEPM * Tmesh^2;
WS2D = WS2D * QEPM * Tmesh^2;

%NOTE: a x-check in the number of photons is here in order, to avoid mistakes in the statistics.
%Given the very low weight (much less than 1 since each time represents an individual photon or electron), 
%the binomial approximation will be generally fine.
%Note that PML will detect slightly more photons through this algorithm (it is interrogated first),
%and will only detect one photon per weighted photon. These approximations are not expected to have
%practical consequences in general.

kL_S1 = 0; kR_S1 = 0; kU_S1 = 0; kD_S1 = 0;
kL_S2 = 0; kR_S2 = 0; kU_S2 = 0; kD_S2 = 0;
tPML_S1(1) = -99999; tPMR_S1(1) = -99999; tPMU_S1(1) = -99999; tPMD_S1(1) = -99999;
tPML_S2(1) = -99999; tPMR_S2(1) = -99999; tPMU_S2(1) = -99999; tPMD_S2(1) = -99999;

for i=1:length(tS1) 
    if     (rand<=WS1L(i)), tPML_S1(kL_S1+1) = tS1(i) + lengthS1L(i)/clight;  kL_S1 = kL_S1+1;  %#ok<AGROW>
    elseif (rand<=WS1R(i)), tPMR_S1(kR_S1+1) = tS1(i) + lengthS1R(i)/clight;  kR_S1 = kR_S1+1;  %#ok<AGROW>
    elseif (rand<=WS1U(i)), tPMU_S1(kU_S1+1) = tS1(i) + lengthS1U(i)/clight;  kU_S1 = kU_S1+1;  %#ok<AGROW>
    elseif (rand<=WS1D(i)), tPMD_S1(kD_S1+1) = tS1(i) + lengthS1D(i)/clight;  kD_S1 = kD_S1+1;  %#ok<AGROW>
    end
end
for i=1:length(tS2)
    if     (rand<=WS2L(i)), tPML_S2(kL_S2+1) = tS2(i) + lengthS2L(i)/clight;  kL_S2 = kL_S2+1;  %#ok<AGROW>
    elseif (rand<=WS2R(i)), tPMR_S2(kR_S2+1) = tS2(i) + lengthS2R(i)/clight;  kR_S2 = kR_S2+1;  %#ok<AGROW>
    elseif (rand<=WS2U(i)), tPMU_S2(kU_S2+1) = tS2(i) + lengthS2U(i)/clight;  kU_S2 = kU_S2+1;  %#ok<AGROW>
    elseif (rand<=WS2D(i)), tPMD_S2(kD_S2+1) = tS2(i) + lengthS2D(i)/clight;  kD_S2 = kD_S2+1;  %#ok<AGROW>
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% GENERATE WAVEFORMS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Generate waveforms (I. Add time entries)
Twvf = -1000:Tbin:10000;
PML_wvf_S1 = zeros(size(Twvf)); PMR_wvf_S1 = zeros(size(Twvf)); PMU_wvf_S1 = zeros(size(Twvf)); PMD_wvf_S1 = zeros(size(Twvf));
PML_wvf_S2 = zeros(size(Twvf)); PMR_wvf_S2 = zeros(size(Twvf)); PMU_wvf_S2 = zeros(size(Twvf)); PMD_wvf_S2 = zeros(size(Twvf));

if(kL_S1>0), PML_wvf_S1 = hist1D(tPML_S1, Twvf); end
if(kR_S1>0), PMR_wvf_S1 = hist1D(tPMR_S1, Twvf); end
if(kU_S1>0), PMU_wvf_S1 = hist1D(tPMU_S1, Twvf); end
if(kD_S1>0), PMD_wvf_S1 = hist1D(tPMD_S1, Twvf); end

if(kL_S2>0), PML_wvf_S2 = hist1D(tPML_S2, Twvf); end
if(kR_S2>0), PMR_wvf_S2 = hist1D(tPMR_S2, Twvf); end
if(kU_S2>0), PMU_wvf_S2 = hist1D(tPMU_S2, Twvf); end
if(kD_S2>0), PMD_wvf_S2 = hist1D(tPMD_S2, Twvf); end

%Generate waveforms (II. Add PM statistics)
for i=1:length(Twvf)
    if(PML_wvf_S1(i)>0),  PML_wvf_S1(i) = random('Normal', PML_wvf_S1(i), sigmaAovA1*sqrt(PML_wvf_S1(i))); end    
    if(PMR_wvf_S1(i)>0),  PMR_wvf_S1(i) = random('Normal', PMR_wvf_S1(i), sigmaAovA1*sqrt(PMR_wvf_S1(i))); end     
    if(PMU_wvf_S1(i)>0),  PMU_wvf_S1(i) = random('Normal', PMU_wvf_S1(i), sigmaAovA1*sqrt(PMU_wvf_S1(i))); end 
    if(PMD_wvf_S1(i)>0),  PMD_wvf_S1(i) = random('Normal', PMD_wvf_S1(i), sigmaAovA1*sqrt(PMD_wvf_S1(i))); end
    
    if(PML_wvf_S2(i)>0),  PML_wvf_S2(i) = random('Normal', PML_wvf_S2(i), sigmaAovA1*sqrt(PML_wvf_S2(i))); end    
    if(PMR_wvf_S2(i)>0),  PMR_wvf_S2(i) = random('Normal', PMR_wvf_S2(i), sigmaAovA1*sqrt(PMR_wvf_S2(i))); end     
    if(PMU_wvf_S2(i)>0),  PMU_wvf_S2(i) = random('Normal', PMU_wvf_S2(i), sigmaAovA1*sqrt(PMU_wvf_S2(i))); end 
    if(PMD_wvf_S2(i)>0),  PMD_wvf_S2(i) = random('Normal', PMD_wvf_S2(i), sigmaAovA1*sqrt(PMD_wvf_S2(i))); end
end
 
%Generate waveforms (III. Convolute S2 with gap response -assumed to be a flat distribution)
responseGap = ones(size(Twvf));
responseGap(Twvf>(Twvf(1)+ gap/vdGap)) = 0;
responseGap = 1/Tbin*responseGap/(sum(responseGap));                        
%NOTE: MTLconvol (and MTLconvolutionTester) not working fine in Matlab2018. It requires explicit normalization!!
PML_wvf_S2 = MTLconvol(PML_wvf_S2, Twvf, responseGap, Twvf, 'tt');
PMR_wvf_S2 = MTLconvol(PMR_wvf_S2, Twvf, responseGap, Twvf, 'tt');
PMU_wvf_S2 = MTLconvol(PMU_wvf_S2, Twvf, responseGap, Twvf, 'tt');
PMD_wvf_S2 = MTLconvol(PMD_wvf_S2, Twvf, responseGap, Twvf, 'tt');

PML_wvf = PML_wvf_S1 + PML_wvf_S2;
PMR_wvf = PMR_wvf_S1 + PMR_wvf_S2;
PMU_wvf = PMU_wvf_S1 + PMU_wvf_S2;
PMD_wvf = PMD_wvf_S1 + PMD_wvf_S2;

%Generate waveforms (IV. Convolute with PM response)
responsePM = exp(-(Twvf - meanTresponse).^2/(2*(sigmaTresponse)));
responsePM = 1/Tbin*responsePM/(sum(responsePM));
PML_wvf = MTLconvol(PML_wvf, Twvf, responsePM, Twvf, 'tt');
PMR_wvf = MTLconvol(PMR_wvf, Twvf, responsePM, Twvf, 'tt');
PMU_wvf = MTLconvol(PMU_wvf, Twvf, responsePM, Twvf, 'tt');
PMD_wvf = MTLconvol(PMD_wvf, Twvf, responsePM, Twvf, 'tt');

%Generate waveforms (V. Add noise, given in photon units)
for i=1:length(Twvf)
    PML_wvf(i) = PML_wvf(i) + random('Normal', 0, sigmaNoise);
    PMR_wvf(i) = PMR_wvf(i) + random('Normal', 0, sigmaNoise);
    PMU_wvf(i) = PMU_wvf(i) + random('Normal', 0, sigmaNoise);
    PMD_wvf(i) = PMD_wvf(i) + random('Normal', 0, sigmaNoise);
end

%Generate waveforms (VI. Convolute with Osci response)
Tau_osci     = 1/BW/2/pi;
responseOsci = 1/Tau_osci*exp(-Twvf/Tau_osci);
responseOsci = 1/Tbin*responseOsci/(sum(responseOsci));
PML_wvf = MTLconvol(PML_wvf, Twvf, responseOsci, Twvf, 'tt');
PMR_wvf = MTLconvol(PMR_wvf, Twvf, responseOsci, Twvf, 'tt');
PMU_wvf = MTLconvol(PMU_wvf, Twvf, responseOsci, Twvf, 'tt');
PMD_wvf = MTLconvol(PMD_wvf, Twvf, responseOsci, Twvf, 'tt');

%Generate waveforms (VII. Convert from photon number to amplitude)
PML_wvf = PML_wvf*A1;
PMR_wvf = PMR_wvf*A1;
PMU_wvf = PMU_wvf*A1;
PMD_wvf = PMD_wvf*A1;

 
