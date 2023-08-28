function dEdx = E_loss_BB(Z, A, rho, z, E, M, delta)

% dEdx   = MeV/cm
% Z      = atomic number of the medium
% A      = mass   number of the medium
% rho    = density [g/cm^3]
% z      = atomic number of the particle
% E      = Energy of the particle [MeV]
% M      = mass of the particle [MeV]
% delta  = Density effect parameter

K   = 0.307075;       %MeV g^-1 cm^2
me  = 0.511;          %MeV
I   = 1e-6*16*Z^0.9;  %Parameterization given in the PDG booklet
%Z   = 9;
%A   = 19;
%z   = 1;
%M  = 0.511;          

gamma   = E/M;
beta    = sqrt(1-1./gamma.^2);

T_max = 2*me*beta.^2.*gamma.^2./(1 + 2*gamma*me/M +(me/M)^2);

dEdx  = K * z^2 * Z/A * 1./beta.^2 .* ( 1./2 * log( 2*me*beta.^2.*gamma.^2.*T_max/I^2 ) - beta.^2);


return

