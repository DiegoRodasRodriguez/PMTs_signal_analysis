%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% script to simulate CMOS sensor reponse (DGD 07/04/19) %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  function [NphXY] = CMOSresponse(xePh, yePh, zePh, tePh, WePh);         %   
%                                                                         %
% INPUT:                                                                  %
% xePh, yePh, zePh, tePh, WePh: position time and weight for S2           %
%                                                                         %
% OUTPUT:                                                                 %
% NphXY: X-Y matrix with nphotons as entries                              %
%                                                                         %
% PARAMETERS:                                                             %
% M         : magnification                                               %
% lensN     : lens number.                                                %
% sigmaNph  : sigma in number of photons per pixel                        %
% QECMOS    : 70% of the 25% of visible light                             %
% rebin     : group pixels according to camera software                   %
% Npixel    : number of pixels per line                                   %
% pixelSize : pixel size [cm]                                             %
%                                                                         %
% TO DO:                                                                  %
% *Include/introduce function for weighted histogram to include WePh      %
% electron by electron                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NphXY = CMOSresponse(xePh, yePh, WePh)

global M lensN sigmaNph QECMOS rebin Npixel pixelSize;
global Tmesh;

%Calculate Solid angle and sensor size
OmegaCMOS = 1/(16*lensN^2) * 1/(1/M + 1)^2;
WCMOS     = Tmesh * OmegaCMOS * QECMOS;

xsideCMOS = Npixel*pixelSize;
ysideCMOS = Npixel*pixelSize;

%Project anode positions
xeCMOS    = M * xePh; %Inverts!
yeCMOS    = M * yePh; %Inverts!

%Create matrix of weights
Xmax = xsideCMOS/2 - rebin*pixelSize; Ymax = ysideCMOS/2 - rebin*pixelSize; %1 bin needs to be subtracted in order to give the right cound
Xmin = -xsideCMOS/2;                  Ymin = -ysideCMOS/2;
NphXY = hist3D([xeCMOS', yeCMOS'], Xmin:rebin*pixelSize:Xmax, Ymin:rebin*pixelSize:Ymax);

sizeNphXY  = size(NphXY); 
NpixXrebin = sizeNphXY(1); 
NpixYrebin = sizeNphXY(2);

%Include statitistical weights and noise
WePhMean = mean(WePh);  %Note: temporary patch
for i=1:NpixXrebin
    for j=1:NpixYrebin
        NphXY(i,j) = random('Poisson', WePhMean*WCMOS*NphXY(i,j)) + random('Normal', 0, sigmaNph*sqrt(rebin*rebin));
    end
end
 
