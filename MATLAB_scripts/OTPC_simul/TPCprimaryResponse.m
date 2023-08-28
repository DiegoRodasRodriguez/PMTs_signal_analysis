%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% script to simulate TPC primary response (DGD 18/03/19) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% function [xe, ye, ze, te, xph, yph, zph, tph, Lph] =                    %
%                               TPCprimaryResponse(x,y,z,t,DeltaE)        %
%                                                                         %
% INPUT:                                                                  %
% x, y, z, t, DeltaE: pos. and time of energy deposit DeltaE [cm, ns, eV] %
%                                                                         %
% OUTPUT:                                                                 %
% xe, ye, ze, te, tph: pos and time of ionization                         %
% tph: time of scintillation [ns]                                         %
% Lph: wavelength of scintillation [nm]                                   %
%                                                                         %
% PARAMETERS                                                              %
% Wi: average energy to create an electron-ion pair [eV]                  %
% Ws: average energy to create a photon             [eV]                  %
% FanoQ: Fano factor for ionization                                       %
% FanoS: Fano factor for scintillation                                    %
% tau1, tau2, Ratio12: time constants (fast+slow) & ratio of F/S          %
%                                                                         %
% NOTES/TO DOs:                                                           %
%                                                                         %
% *Include scintillation spectrum                                         %
% *Fano included for constant energy loss (modify when B-B included!)     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [xe, ye, ze, te, xph, yph, zph, tph, Lph]  = TPCprimaryResponse(x, y, z, t, DeltaE)

global Wi Ws tau1 tau2 Ratio12 FanoQ FanoS;

neToTmean = round(sum(DeltaE)/Wi); nphToTmean = round(sum(DeltaE)/Ws);

sigmaStepQ = ( sqrt(FanoQ) * sqrt(neToTmean)  ) / sqrt(length(DeltaE));
sigmaStepS = ( sqrt(FanoS) * sqrt(nphToTmean) ) / sqrt(length(DeltaE));

ne    = zeros(1, length(DeltaE));
nph   = zeros(1, length(DeltaE));

neToT = 0; nphToT = 0;

for i=1:length(DeltaE)    
    ne_mean   = DeltaE(i)/Wi;
    ne(i)     = round(random('Normal', ne_mean,  sigmaStepQ));     
    nph_mean  = DeltaE(i)/Ws;
    nph(i)    = round(random('Normal', nph_mean, sigmaStepS)); 
    neToT     = neToT  + ne(i); 
    nphToT    = nphToT + nph(i); 
end

xe  = zeros(1,neToT);  ye  = zeros(1,neToT);  ze  = zeros(1,neToT);  te  = zeros(1,neToT);
xph = zeros(1,nphToT); yph = zeros(1,nphToT); zph = zeros(1,nphToT); tph = zeros(1,nphToT);
Lph = zeros(1,nphToT);

I0 = 1; nphCumm = 0;

%Loop over energy deposits
for i=1:length(DeltaE)    
   %Electrons
   If        = I0 + ne(i);
   xe(I0:If) = x(i);
   ye(I0:If) = y(i);
   ze(I0:If) = z(i);
   te(I0:If) = t(i);
   I0        = If + 1;
   
   %Photons
   for j=1:nph(i)
       if(rand>Ratio12), tph(nphCumm+j) = random('Exponential', tau2);
       else,             tph(nphCumm+j) = random('Exponential', tau1);
       end
       xph(nphCumm+j) = x(i);
       yph(nphCumm+j) = y(i);
       zph(nphCumm+j) = z(i);
       Lph(nphCumm+j) = -1;
   end
   nphCumm   = nphCumm + nph(i);
end
