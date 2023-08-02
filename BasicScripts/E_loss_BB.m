function dEdx = E_loss_BB(Z_A, I, rho, z, Ekin, M, delta_o, C, X_o, X_1, a, m)

% dEdx   = keV/cm

% Z      = atomic number of the medium
% A      = mass   number of the medium
% rho    = density [g/cm^3]
% z      = atomic number of the particle
% E      = Energy of the particle [MeV]
% M      = mass of the particle [MeV]
% delta  = Density effect parameter

K   = 307.075;        %keV g^-1 cm^2
me  = 0.511;          %MeV
I   = I*1e-6;         %MeV
%Z   = 9;
%A   = 19;
%z   = 1;
%M  = 0.511;          

E = Ekin + M; 

gamma   = E./M;
beta    = sqrt(1-1./gamma.^2);
gbeta   = gamma.*beta;

T_max   = 2*me*beta.^2.*gamma.^2./(1 + 2*gamma*me/M +(me/M)^2);

%density effect. From Sternheimer et al. Atomic and Nuclear Data Tables
%30(1984)261

I1 = find(log10(gbeta)<=X_o);
I2 = find(log10(gbeta)>X_o & log10(gbeta)<X_1);
I3 = find(log10(gbeta)>=X_1);

gbeta1 = gbeta(I1);
gbeta2 = gbeta(I2);
gbeta3 = gbeta(I3);

delta1 = delta_o*10.^(2*(log10(gbeta1)-X_o));
delta2 = 4.6052*log10(gbeta2) + a*(X_1-log10(gbeta2)).^m - C;
delta3 = 4.6052*log10(gbeta3) - C;

delta = [delta1 delta2 delta3];

%computing energy loss

dEdx  = rho*K * z^2 * Z_A * 1./beta.^2 .* ( 1./2 * log( 2*me*beta.^2.*gamma.^2.*T_max/I^2 ) - beta.^2 - delta/2);

return

