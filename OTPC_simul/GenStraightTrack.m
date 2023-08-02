%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% script to simulate straight tracks (DGD 18/03/19) %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
% function [x,y,z,t,DeltaE] =                                             %
%            GenStraightTrack(z0,theta0,phi0,Range,dEdx,M,nsteps)         %
%                                                                         %
% INPUT:                                                                  %
% z0    : initial z position [cm]                                         %
% theta0: initial theta      [deg]                                        %
% phi0  : initial phi        [deg]                                        %
% Range : track length       [cm]                                         %
% dEdx  : energy loss        [MeV/cm]                                     %
% M     : particle mass      [MeV]                                        %
% nsteps: n of energy deposits                                            %
%                                                                         %
% OUTPUT:                                                                 %
% x, y, z, t, DeltaE: pos. and time of energy deposit DeltaE [cm, ns, eV] %
%                                                                         %
% NOTES/TO DOs:                                                           %
% *Include Bethe-Bloch                                                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,y,z,t,DeltaE] = GenStraightTrack(z0, theta0, phi0, Range, dEdx, M, nsteps)

Ekin    = dEdx*Range;          % [MeV]
clight  = 30;                  % [cm/ns]

rlength = Range*rand(1,nsteps);
rlength = sort(rlength);
beta    = sqrt(2*Ekin/M);
vel     = beta*clight;

t         = zeros(size(rlength));
totLength = 0;
for i=2:length(rlength)
    deltaLength = rlength(i) - rlength(i-1);    
    totLength   = totLength + deltaLength;
    t(i)        = t(i-1) + deltaLength/vel;
    beta        = sqrt(2*(Ekin - dEdx*totLength)/M);
    vel         = beta*clight;    
end

x       = rlength * sin(theta0*pi/180) * cos(phi0*pi/180);
y       = rlength * sin(theta0*pi/180) * sin(phi0*pi/180);
z       = z0 + rlength * cos(theta0*pi/180);
DeltaE  = Ekin/nsteps*ones(size(rlength)) * 1e+6 ;

end
