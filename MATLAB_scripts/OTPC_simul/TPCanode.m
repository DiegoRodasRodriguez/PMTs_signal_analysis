%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% script to simulate anode response (DGD 18/03/19) %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%  function [xePh, yePh, zePh, tePh, WPh] = TPCanode(xeA, yeA, zeA, teA)  %
%                                                                         %
% INPUT:                                                                  %
% xeA, yeA, zeA, teA: pos. and time of electrons at anode                 %
%                                                                         %
% OUTPUT:                                                                 %
% xePh, yePh, zePh, tePh: pos and time of photons produced at anode       %
% WPh:  statistical weight of photons produced at anode (only one         %
%       'effective' photon is created per electron)                       %
%                                                                         %
% PARAMETERS:                                                             %
% hole    : hole diameter      [cm]                                       %
% pitch   : pitch size         [cm]                                       %
% gap     : thickness of scintillation region [cm]                        %
% OptGain : optical gain                                                  %
% vdGap   : velocity across gap [cm/ns]                                   %
% tau1Gap, tau2Gap, Ratio12Gap: time constants (fast+slow) & ratio of F/S %
%                                                                         %
% NOTES/TO DOs:                                                           %
% *It assumes Furry law to calculate statistical weight, assigned to anode%
% *Scintillation during electron transit included through a convolution in%
% PMTresponse (large numerical overhead otherwise)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xePh, yePh, zePh, tePh, WePh] = TPCanode(xeA, yeA, zeA, teA)

global hole pitch gap OptGain vdGap tau1Gap tau2Gap Ratio12Gap; %#ok<NUSED>

xePh = pitch*round(xeA/pitch);
yePh = pitch*round(yeA/pitch);
zePh = zeA;

tePh = zeros(size(zeA)); WePh = zePh;

for i=1:length(xeA)
    %smear within the hole
    r2_rnd  = rand;    
    r_rnd   = sqrt(r2_rnd)*hole/2;
    phi_rnd = 2*pi*rand;
    xePh(i) = xePh(i) + r_rnd*sin(phi_rnd);    
    yePh(i) = yePh(i) + r_rnd*cos(phi_rnd);
    %Add avalanche statistics
    WePh(i) = random('exponential', OptGain);   
    if(rand>Ratio12Gap), tePh(i) = teA(i) + random('Exponential', tau2Gap);
    else,                tePh(i) = teA(i) + random('Exponential', tau1Gap);
    end
end

