function MinFun = GFitterMERT2swarm(par)

global Ydata    Xdata    STDdata    ERRdata    STAdata
global Ydata_DT Xdata_DT STDdata_DT ERRdata_DT STAdata_DT
global Ydata_DL Xdata_DL STDdata_DL ERRdata_DL STAdata_DL

global OptFit Nseries isQinFit isErrorType 
global Qm_KuroFit Qme_KuroFit ener_KuroFit
global Qm0 e0 edNde
global e_match lambda
global E_scan EN_scan N M P T Z isthermOff DLverbose

Qleft    = MERT_Qm_for_fit(par, edNde);
Qright   = Qm0;
Qright   = interp1(e0, Qm0, edNde);
Qm       = Qmatcher(Qleft, Qright, edNde, e_match, lambda);

Qmodel   = MERT_Qm_for_fit(par, ener_KuroFit);

Nmu_model      = zeros(size(E_scan));
Nmu_data       = zeros(size(E_scan));
Nmu_dataERROR  = ones(size(E_scan));

NDT_model      = zeros(size(E_scan));
NDT_data       = zeros(size(E_scan));
NDT_dataERROR  = ones(size(E_scan));

NDL_model      = zeros(size(E_scan));
NDL_data       = zeros(size(E_scan));
NDL_dataERROR  = ones(size(E_scan));

if    (OptFit == 0), isVdFit = 1; isDTFit = 0; isDLFit = 0; isQFit = 0; Nseries = 1;
elseif(OptFit == 1), isVdFit = 1; isDTFit = 1; isDLFit = 0; isQFit = 0; Nseries = 2;
elseif(OptFit == 2), isVdFit = 0; isDTFit = 1; isDLFit = 0; isQFit = 0; Nseries = 1;
elseif(OptFit == 3), isVdFit = 0; isDTFit = 0; isDLFit = 1; isQFit = 0; Nseries = 1;
elseif(OptFit == 3), isVdFit = 1; isDTFit = 0; isDLFit = 1; isQFit = 0; Nseries = 2;
elseif(OptFit == 5), isVdFit = 0; isDTFit = 1; isDLFit = 1; isQFit = 0; Nseries = 2;
elseif(OptFit == 6), isVdFit = 1; isDTFit = 1; isDLFit = 1; isQFit = 0; Nseries = 3;
elseif(OptFit == 7), isVdFit = 1; isDTFit = 0; isDLFit = 0; isQFit = 1; Nseries = 1;
elseif(OptFit == 8), isVdFit = 1; isDTFit = 1; isDLFit = 0; isQFit = 1; Nseries = 2;
elseif(OptFit == 9), isVdFit = 1; isDTFit = 1; isDLFit = 1; isQFit = 1; Nseries = 3;
end

%NOTE: if error == 0 complete with statistical estimate (needs only one measurement)

for i=1:length(E_scan)    
   if(isVdFit)
        vd              = vdCALC_optim(Qm, edNde, edNde, M, E_scan(i), P, T, Z, isthermOff);
        Nmu_model(i)    = vd * N / (E_scan(i) * 100) / 100 / 1e+22;
        Nmu_data(i)     = Ydata(min(find(Xdata>=EN_scan(i))));
        if    (isErrorType==1), Nmu_dataERROR(i) = STAdata(min(find(Xdata>=EN_scan(i))));
        elseif(isErrorType==2), Nmu_dataERROR(i) = ERRdata(min(find(Xdata>=EN_scan(i))));
        elseif(isErrorType==3), Nmu_dataERROR(i) = STDdata(min(find(Xdata>=EN_scan(i))));
        end
        if(Nmu_dataERROR(i) == 0), Nmu_dataERROR(i) = STAdata(min(find(Xdata>=EN_scan(i)))); end
   end
   if(isDTFit)
       DT               = DTCALC_optim(Qm, edNde, edNde, M, E_scan(i), P, T, Z, isthermOff);
       NDT_model(i)     = DT * N / 1e+22 * 1e-2;
       NDT_data(i)      = Ydata_DT(min(find(Xdata>=EN_scan(i))));
       if    (isErrorType==1), NDT_dataERROR(i) = STAdata_DT(min(find(Xdata_DT>=EN_scan(i))));
       elseif(isErrorType==2), NDT_dataERROR(i) = ERRdata_DT(min(find(Xdata_DT>=EN_scan(i))));
       elseif(isErrorType==3), NDT_dataERROR(i) = STDdata_DT(min(find(Xdata_DT>=EN_scan(i))));
       end
       if(NDT_dataERROR(i) == 0), NDT_dataERROR(i) = STAdata_DT(min(find(Xdata_DT>=EN_scan(i)))); end
   end
   if(isDLFit)
       DL               = DLCALC_optim(Qm, edNde, edNde, M, E_scan(i), P, T, Z, isthermOff, vd, DT, DLverbose);
       NDL_model(i)     = DL * N / 1e+22 * 1e-2;
       NDL_data(i)      = Ydata_DL(min(find(Xdata>=EN_scan(i))));
       if    (isErrorType==1), NDL_dataERROR(i) = STAdata_DL(min(find(Xdata_DL>=EN_scan(i))));
       elseif(isErrorType==2), NDL_dataERROR(i) = ERRdata_DL(min(find(Xdata_DL>=EN_scan(i))));
       elseif(isErrorType==3), NDL_dataERROR(i) = STDdata_DL(min(find(Xdata_DL>=EN_scan(i))));
       end       
       if(NDL_dataERROR(i) == 0), NDL_dataERROR(i) = STAdata_DL(min(find(Xdata_DL>=EN_scan(i)))); end
   end
end


MinFun  = sqrt(isVdFit*((Nmu_model - Nmu_data).^2)./(Nmu_dataERROR.^2) + isDTFit*((NDT_model - NDT_data).^2)./(NDT_dataERROR.^2) + isDLFit*((NDL_model - NDL_data).^2)./(NDL_dataERROR.^2) + ...
    sum(isQinFit*((Qmodel - Qm_KuroFit).^2)./(Qme_KuroFit.^2))/length(Qmodel)*ones(size(Nmu_model))); %Calculates the total Chi2 re-weighted to the same number of points as the others  

% disp(['MinFun = ', num2str(sum(MinFun))]);
% disp(['Qsum = ', num2str(sum(sum(isQinFit*((Qmodel - Qm_KuroFit).^2)./(Qme_KuroFit.^2))/length(Qmodel)*ones(size(Nmu_model))))]);
end
