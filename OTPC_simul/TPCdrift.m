%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% script to simulate TPC drift (DGD 18/03/19) %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%       function [xeA, yeA, zeA, teA] = TPCdrift(xe, ye, ze, te)          %
%                                                                         %
% INPUT:                                                                  %
% xe, ye, ze, te: pos. and time of initial electrons                      %
%                                                                         %
% OUTPUT:                                                                 %
% xeA, yeA, zeA, teA: pos and time of electrons at anode                  %
%                                                                         %
% PARAMETERS:                                                             %
% vd    : drift velocities       [cm/ns]                                  %
% DL, DT: diffusion coefficient  [cm/sqrt(cm)]                            %
% P     : operating P            [bar]                                    %
% eta   : attachment coefficient [cm^-1]                                  %
%                                                                         %
% NOTES/TO DOs:                                                           %
% Include initiatialization parameters from file                          %
% Add attachment                                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xeA, yeA, zeA, teA] = TPCdrift(xe, ye, ze, te)

global P vd DL DT; %#ok<NUSED>

xeA = zeros(size(xe)); yeA = xeA; zeA = xeA; teA = xeA;

for i=1:length(xe)
    teA(i) = te(i) + random('Normal', ze(i)/vd, (DL*sqrt(ze(i)))/vd);
    xeA(i) = random('Normal', xe(i),    DT*sqrt(ze(i)));
    yeA(i) = random('Normal', ye(i),    DT*sqrt(ze(i)));
    zeA(i) = 0;
end
 
