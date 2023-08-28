%FIXME: include rate drop along perpendicular direction
%NOTE: Analysis for obtaining the profiles must be re-done in order
%to avoid this nasty shift.

% HV for strip position scan: 6.45 kV (JB)          -> 103 kV/cm (too close to the plateau drop)
% HV for pad rate scan:       ~5.6 kV (from Figure) -> 101 kV/cm (also well on plateau)
% HV for pad0 rate scan:      ~6.6 kV (paper)       -> 120 kV/cm. Well on the plateau

clear all;
close all;

%FIRST PLOT EFFICIENCIES AND Q vs E relations

%Charge vs E at low rates (from .xls)
Q_strip   = [1.45 2.15 2.98  3.80  4.72 5.74  6.67  7.96 9.50 11.67]/2;
E_strip   = [5.6   5.8  6.0  6.2   6.4  6.6   6.7  6.9   7    7.1]*2/5/0.25*10;
eff_strip = [0.25 0.64 0.875 0.95 0.975 0.975 0.98 0.99 0.99 0.99];
%Charge vs E at low rates (scanned from figure)
Q_pad    = [1.5 2 2.5 3.15 3.5 3.75 4.4 4.8 6.5]/2;
E_pad    = [4.8    5   5.2  5.4  5.6  5.8  5.9  6.0 6.2]*2/5/0.22*10;
eff_pad  = [0.76 0.83 0.85 0.86 0.86 0.86 0.86 0.86 0.86]/0.86;
%eff(up-cut) vs E for strip (from .xls)
eff_pad0 = [0.8 0.865 0.9 0.925 0.935 0.95 0.95 0.955 0.96 0.955 0.96 0.955 0.96]/0.96;
E_pad0   = [5.8 5.9   6.0 6.1   6.2   6.3  6.4  6.5   6.6  6.7   6.8  6.9   7.0]*2/5/0.22*10;

%
E_ana         = 85:0.01:130;
alfa          = GASalphaFreon(E_ana)*0.85;
eta           = 0;
n_o_mips      = 7.5;
Npad09        = n_o_mips*5*2*0.22;
Npad11        = 1.5*n_o_mips*5*2*0.22;
Nstrip11      = 1.5*n_o_mips*5*2*0.25;

isSignalDiff  = 0;                          %NOTE: at the moment only single-ended geometries are studied. This could become important.
Zdet          = 15;                         %FIXME.
tau_strip     = 2*Zdet/(Zdet+50);
Eweight_ideal = 0.53;
Qth           = 0.055;
%%%%%%%%%%%%%%%%%%%%%%
eff_ana_Pad_GSI09   = PS2EffAnalytic(alfa, eta, Npad09,   0.22, Qth, Eweight_ideal, 1, 0);
eff_ana_Pad_GSI11   = PS2EffAnalytic(alfa/0.85*1.20, eta, Npad11,   0.22, Qth, Eweight_ideal, 1, 0);
eff_ana_Strip_GSI11 = PS2EffAnalytic(alfa/0.85*1.25, eta, Nstrip11, 0.25, Qth, Eweight_ideal*tau_strip, 0, 0);

%Efficiency

figure;
plot(E_strip, eff_strip, 'o'); hold on;
plot(E_pad,   eff_pad,   'ro');
plot(E_pad0,  eff_pad0,  'go');
plot(E_ana, eff_ana_Pad_GSI09,   'g-');
plot(E_ana, eff_ana_Pad_GSI11,   'r-');
legend('strip HZDR11', 'pad HZDR11', 'pad GSI09'); yaxis(0,1);
plot(E_ana, eff_ana_Strip_GSI11, '-'); hold on;
xlabel('E[kV/cm]'); ylabel('efficiency');
%%%%%%%%%%%%%%%%%%%
%% Q vs E relations
%%%%%%%%%%%%%%%%%%%

figure;
plot(E_strip, Q_strip, 'o'); hold on;
plot(E_pad,   Q_pad,   'sr'); 
xlabel('E [kV/cm]');
ylabel('average charge per gap (I/r/2) [pC]');
legend('strip', 'pad');

p_strip = polyfit(E_strip(find(E_strip<107)), Q_strip(find(E_strip<107)), 1);
plot(E_strip(E_strip<107), polyval(p_strip,E_strip(E_strip<107)), '-b');
a_strip   =  p_strip(1);             % [pC/(kV/cm)]
Eth_strip = -p_strip(2)/p_strip(1);  % [kV/cm]

p_pad   = polyfit(E_pad(find(E_pad<110)), Q_pad(find(E_pad<110)), 1);
plot(E_pad(E_pad<110), polyval(p_pad,E_pad(E_pad<110)), '-r');
a_pad   =  p_pad(1);            % [pC/(kV/cm)]
Eth_pad = -p_pad(2)/p_pad(1);   % [kV/cm]

%DC-model
flux      = 0.01:0.01:300;        % [kHz/cm^2]
AreaIrrad = 3*3;                  % [cm^2]
AreaEff   = 3*3;                  % [cm^2]
dglass_average   = 6*0.07/5;      % [cm]
rho_plate        = 2e+10;         % [Ohms x cm]

Reff             = rho_plate*dglass_average/AreaEff; 
E_strip_         = 6.45*2/5/0.25*10;

[Ihv, Vgap, Egap, qgap] = PS2_lineal_DCmodel(flux, AreaIrrad, E_strip_, a_strip, Eth_strip, Reff, 0.25, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rate, Flux and position scan
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rateS1 = [58 160 230 410 578 826 979 1348 1602 1890];

load F:\HOME\_PEOPLE\_JINGBO\HZDR_tests_April\Profile_vs_position
X=U(:,1)*2.31; Flux_=U(:,2:11);  Flux_(1,:)=0;
Xoffset =  10.8;
Xshift  = X - Xoffset;
for i=1:10
    I   = find(Xshift>0);
    XX  = Xshift(I);
    FF  = Flux_(:,i);
    FF  = FF(I);
    FF_ = interp1(XX, FF, X(1:length(XX)-1));
    FF_(length(X)) = 0; 
    Flux(:,i) = FF_;
end
Err_Flux = sqrt(Flux);

load F:\HOME\_PEOPLE\_JINGBO\HZDR_tests_April\ProfileU_vs_position
X=U(:,1)*2.31; Flux_U=U(:,2:11); Err_Flux_U = sqrt(Flux_U); Flux_U(1,:)=0;

for i=1:10
    I   = find(Xshift>0);
    XX  = Xshift(I);
    FF  = Flux_U(:,i);
    FF  = FF(I);
    FF_ = interp1(XX, FF, X(1:length(XX)-1));
    FF_(length(X)) = 0; 
    Flux_U(:,i) = FF_;
end

load F:\HOME\_PEOPLE\_JINGBO\HZDR_tests_April\Qtot_vs_position
X=U(:,1)*2.31; Qtot=U(:,2:11)/2/1000;

load F:\HOME\_PEOPLE\_JINGBO\HZDR_tests_April\V_vs_position
X=U(:,1)*2.31; HV=U(:,2:11);

load F:\HOME\_PEOPLE\_JINGBO\HZDR_tests_April\eff_vs_position
X=U(:,1)*2.31; Eff=U(:,2:11)/100;

%SOME PARAMETERS USED LATER ON
Xmin = 80; Xmax = 160;
Icut = find(X>Xmin & X<Xmax);
XHistbin    = 0.25; %In cm
XDistbin    = 0.01; %In cm
strip_width = 2.2;
pitch       = 2.5;
Area        = XHistbin*strip_width;  

%LINEAR PROFILES
NormXRPChist    = 1/(sum(Flux(:,1))*XHistbin);          %Integral
NormXRPChist_U  = 1/(sum(Flux_U(:,1))*XHistbin);
NormXScindist   = 1/sqrt(2*pi*1.1^2);
NormYScindist   = 1/sqrt(2*pi*1.5^2);
NormXRPCdist    = 0.17;                 %From fit. Not normalized to one, strictly.
NormYRPCdist    = 1/sqrt(pi*3.186^2);

Xgg = -12:XDistbin:12;
Ygg = -12:XDistbin:12;

%Histograms    [cm^-1]
dPdXhist_RPCFlux    = NormXRPChist*Flux(:,1);
dPdXhist_RPCFlux_U  = NormXRPChist_U*Flux_U(:,1);
dP_Err_Flux         = NormXRPChist_U*Err_Flux_U(:,1);
%Functional distributions [cm^-1]
dPdX_Scin           = NormXScindist*exp(-(Xgg/sqrt(2)/1.1).^2);
dPdY_Scin           = NormYScindist*exp(-(Ygg/sqrt(2)/1.5).^2);
dPdX_RPCself        = NormXRPCdist *exp(-(Xgg/2.81).^2);
dPdY_RPCself_extrap = NormYRPCdist *exp(-(Xgg/3.186).^2);
%Other functional distributions
dPdX_wRPC    = 1/sqrt(2*pi*(1.1^2+0.5^2))*exp(-(Xgg/sqrt(2)/sqrt(1.1^2+0.5^2)).^2);
dPdX_RPCtrig = 0.25*exp(-(Xgg/2.041).^2);

%%%LINEAR PROBABILITY DISTRIBUTIONS [cm^-1]
figure; 
subplot(1,2,1);%(along X)
XcentFit = 124.1;
plot((X - XcentFit + Xoffset)/10, dPdXhist_RPCFlux, 's'); hold on;
plot((X - XcentFit + Xoffset)/10, dPdXhist_RPCFlux_U, 'or');
plot(Xgg, dPdX_Scin,       '-g');
plot(Xgg, dPdX_wRPC, '--g');
legend('along 2.2cm-wide strip (external trigger)', 'along 2.2cm-wide strip (self-trigger)', 'from finger scintillator over the RPC box (self-trigger)', 'from finger scintillator after including the RPC position resolution');
plot(Xgg, dPdX_RPCself, '-r');
plot(Xgg, dPdX_RPCtrig, '-');
errorbar((X-XcentFit + Xoffset)/10, dPdXhist_RPCFlux_U, dP_Err_Flux, 'or');
xaxis(-6,6); yaxis(0,0.5); xlabel('position[cm]'); ylabel('dP/dX [cm^{-1}]');
subplot(1,2,2);%(along Y)
plot(Ygg, dPdY_RPCself_extrap, '-r'); hold on;
plot(Ygg, dPdY_Scin, '-g');
legend('frpm the RPC (extrapolated)', 'from finger scintillator over the RPC box (self-trigger)');

%If one assumes that X-shape and Y-shape are independent from Y and X,
%respectively. One can take any of the profiles as representative of the
%integrated profile. So, it is possible to calculate the linear flux
%(integrated over the perpendicular direction).

%%%LINEAR FLUX [kHz/cm]
figure;
subplot(1,2,1);
plot((X-XcentFit + Xoffset)/10, rateS1(1)*dPdXhist_RPCFlux, 's'); hold on;
plot((X-XcentFit + Xoffset)/10, rateS1(1)*dPdXhist_RPCFlux_U, 'or');
plot(Xgg, rateS1(1)*dPdX_Scin, '-g');
plot(Xgg, rateS1(1)*dPdX_wRPC, '--g');
legend('RPC ext trigg', 'RPC self trigg', 'scint self trigg', 'scint self trigg + RPC res');
title('parallel to floor [x] (S_{21} rate =58kHz)');
xaxis(-6,6); yaxis(0,(rateS1(1)/sqrt(2*pi*1.1^2))*1.2); xlabel('position[cm]'); ylabel('\Phi_L (integrated along y) [kHz/cm]');
subplot(1,2,2)
plot(Ygg, rateS1(1)*dPdY_RPCself_extrap, '-r'); hold on;
plot(Ygg, rateS1(1)*dPdY_Scin, '-g');
legend('RPC self trigg (extrp)', 'scint self trigg');
title('perpendicular to floor [x] (S_{21} rate =58kHz)');
xaxis(-6,6); yaxis(0,(rateS1(1)/sqrt(2*pi*1.5^2))*1.2); xlabel('position[cm]'); ylabel('\Phi_L (integrated along x) [kHz/cm]');

%%%2D FLUX [kHz/cm]
figure;

Xs=[-1.1 -1.1 1.1 1.1]; Ys=[-12 12 12 -12]; Zs=[0 0 0 0];
title('Distribution from scintillator');
plot3(Xgg, zeros(size(Xgg)), rateS1(1)*NormXScindist*dPdY_Scin, 'k'); hold on;
plot3(zeros(size(Ygg)), Ygg, rateS1(1)*NormYScindist*dPdX_Scin, 'k');

fill3(Xs, Ys, Zs, 'b');
fill3(Xs-2.5, Ys, Zs, 'b');
fill3(Xs+2.5, Ys, Zs, 'b');
box; zlabel('Flux [kHz/cm^2]'); xaxis(-4,4); yaxis(-12,12);

figure;

Xs=[-1.1 -1.1 1.1 1.1]; Ys=[-12 12 12 -12]; Zs=[0 0 0 0];
title('Distribution from RPC');
plot3(Xgg, zeros(size(Xgg)), rateS1(1)*NormXRPCdist*dPdY_RPCself_extrap, 'k'); hold on;
plot3(zeros(size(Ygg)), Ygg, rateS1(1)*NormYRPCdist*dPdX_RPCself,        'k');

fill3(Xs, Ys, Zs, 'b');
fill3(Xs-2.5, Ys, Zs, 'b');
fill3(Xs+2.5, Ys, Zs, 'b');
box; zlabel('Flux [kHz/cm^2]'); xaxis(-4,4); yaxis(-12,12);

%Normalization along the Y axis (average over strip)
NormScin2D = rateS1(1)*(sum(dPdY_Scin(find(abs(Ygg)<pitch/2))))*XDistbin/pitch;
NormRPC2D  = rateS1(1)*(sum(dPdY_RPCself_extrap(find(abs(Ygg)<pitch/2))))*XDistbin/pitch;

Flux_scinY_ref_     = NormScin2D*dPdX_Scin;
Flux_scinY_ref(:,1) = interp1(Xgg, Flux_scinY_ref_, (X + Xoffset - XcentFit)/10);

Flux_RPCY_ref_      = NormRPC2D*dPdX_RPCself;
Flux_RPCY_ref(:,1)  = interp1(Xgg, Flux_RPCY_ref_, (X + Xoffset - XcentFit)/10);

figure;
plot(X, Flux_scinY_ref(:,1), '-');      hold on;
plot(X, Flux_RPCY_ref(:,1), '-r');
ylabel('Flux[kHz/cm^2]'); xlabel('position[mm]'); legend('Flux from scintillator', 'Flux from RPC');

for i=1:10
    %Area of the central strip (here the condition requested is that the charge of strip 1 is higher [Check with JB] (probably not done like that.FIXME).
    if(i>1) 
        Flux_scinY_ref(:,i) = Flux_scinY_ref(:,1)*rateS1(i)/rateS1(1); 
        Flux_RPCY_ref(:,i)  = Flux_RPCY_ref(:,1) *rateS1(i)/rateS1(1); 
    end
    
    Norm           = sum(Flux(:,i));
    Flux_Norm(:,i) = Flux(:,i)/Norm * rateS1(i)/Area;
    Err_Flux(:,i)  = Err_Flux(:,i)/Norm * rateS1(i)/Area;
    Flux_Ref(:,i)  = Flux_Norm(:,1)*rateS1(i)/rateS1(1);

    Norm_U           = sum(Flux_U(:,i));
    Flux_U_Norm(:,i) = Flux_U(:,i)/Norm_U * rateS1(i)/Area;    
    Err_Flux_U(:,i)  = Err_Flux_U(:,i)/Norm_U * rateS1(i)/Area;
    Flux_U_Ref(:,i) = Flux_U_Norm(:,1)*rateS1(i)/rateS1(1);
    
    Icut = find(X>80 & X<160);
    Flux_Norm_ = Flux_Norm(:,i); Flux_Ref_ = Flux_Ref(:,i); 
    Flux_scinY_ref_ = Flux_scinY_ref(:,i); Flux_RPCY_ref_ = Flux_RPCY_ref(:,i);
    Flux_U_Ref_ = Flux_U_Ref(:,i); Qtot_ = Qtot(:,i); Eff_ = Eff(:,i);

    Flux_Norm_Cut(:,i)  = Flux_Norm_(Icut);
    Flux_Ref_Cut(:,i)   = Flux_Ref_(Icut);
    Flux_U_Ref_Cut(:,i) = Flux_U_Ref_(Icut);
    Flux_scinY_ref_Cut(:,i) = Flux_scinY_ref_(Icut);    
    Flux_RPCY_ref_Cut(:,i)  = Flux_RPCY_ref_(Icut);
    Qtot_Cut(:,i)       = Qtot_(Icut);
    Eff_Cut(:,i)        = Eff_(Icut);
    X_Cut               = X(Icut);
end

figure; plot(X, Flux_Norm, '-');   ylabel('RPC Flux[kHz/cm^2]');                                            xlabel('position[mm]'); xaxis(50,220);
figure; plot(X, Flux_U_Norm, '-'); ylabel('RPC Flux[kHz/cm^2]');                                            xlabel('position[mm]'); xaxis(50,220);

figure; plot(X, Flux_scinY_ref, '-'); ylabel('scintillator Flux (averaged over one pitch) [kHz/cm^2]');     xlabel('position[mm]'); xaxis(50,220);

figure; plot(X, Flux_Ref, '-');    ylabel('RPC Flux (extrapolated from low intensity) [kHz/cm^2]');         xlabel('position[mm]'); xaxis(50,220);
figure; plot(X, Qtot, '-');        ylabel('Average charge per gap [pC]');                                   xlabel('position[mm]'); xaxis(50,220); yaxis(1,3);
figure; plot(X, Eff , '-');        ylabel('extrapolated efficiency');                                       xlabel('position[mm]'); xaxis(50,220); yaxis(0.7,1.1);

figure; plot(X_Cut, Flux_Norm_Cut, '-'); ylabel('RPC Flux[kHz/cm^2]');                                      xlabel('position[mm]'); xaxis(Xmin,Xmax);
figure; plot(X_Cut, Flux_Ref_Cut, '-');  ylabel('RPC Flux (extrapolated from low intensity) [kHz/cm^2]');   xlabel('position[mm]'); xaxis(Xmin,Xmax);
figure; plot(X_Cut, Qtot_Cut, '-');      ylabel('Average charge per gap [pC]');                             xlabel('position[mm]'); xaxis(Xmin,Xmax); yaxis(1,3);
figure; plot(X_Cut, Eff_Cut , '-');      ylabel('extrapolated efficiency');                                 xlabel('position[mm]'); xaxis(Xmin,Xmax); yaxis(0.7,1.1);

figure;
plot(Qtot_Cut, Eff_Cut, 'o'); xlabel('total charge per gap [pC]'); ylabel('efficiency');

figure;
plot(Flux_Ref_Cut, Qtot_Cut, 'o'); logx; hold on;
%DC-model
flux      = 0.01:0.01:300;        % [kHz/cm^2]
AreaIrrad = 3*3;                  % [cm^2]
AreaEff   = 3*3;                  % [cm^2]
dglass_average   = 6*0.07/5;      % [cm]
rho_plate        = 3.5e+10;       % [Ohms x cm]

Reff      = rho_plate*dglass_average/AreaEff; 
E_strip   = 6.45*2/5/0.25*10;

[Ihv, Vgap, Egap, qgap] = PS2_lineal_DCmodel(flux, AreaIrrad, E_strip, a_strip, Eth_strip, Reff, 0.25, 1);

plot(flux, qgap*0.95, '-'); xlabel('Flux [kHz/cm^2]'); ylabel('total charge [pC]'); %FIXME: this 0.95 must be changed.
figure;
plot(Flux_scinY_ref_Cut, Qtot_Cut, 'o'); logx; hold on;
plot(flux, qgap*0.95, '-'); xlabel('Flux [kHz/cm^2]'); ylabel('total charge [pC]'); %FIXME: this 0.95 must be changed.

figure;
plot(Flux_RPCY_ref_Cut, Qtot_Cut, 'o'); logx; hold on;
plot(flux, qgap*0.95, '-'); xlabel('Flux [kHz/cm^2]'); ylabel('total charge [pC]'); %FIXME: this 0.95 must be changed.

% figure;
% plot(Flux_scinY_ref_Cut(:,1), Qtot_Cut(:,1), 'o'); logx; hold on;

% 
% Flux_scinY_ref_Cut
% 
% plot(Flux_U_Ref_Cut, Qtot_Cut, 'o'); logx; hold on;



%plot(flux, qgap*0.95, '-'); xlabel('Flux [kHz/cm^2]'); ylabel('total charge [pC]'); %FIXME: this 0.95 must be changed.


mosaic; return;