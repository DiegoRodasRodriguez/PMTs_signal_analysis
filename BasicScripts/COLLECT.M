load collect

if exist('Sigmas') % called after anabig.m
   d=[str2num(filename(5:8)) Sigmas' SigmaT0T1 Eff_TOF_Hardware Eff_after_external_cuts MeanA ...
      MeanT1 MeanDTrpc/1000 MeanT2 SlopeTA2 TailsRight TailsLeft WalkT1 WalkT2 WalkArpc SigmaRpcT1uncorr];
end

if exist('N12NoGood') % called after Yeff.m
   d=[str2num(filename(5:8)) Events N1Good N2Good N12Good N12NoGood TimeEff1 TimeEff2];   
end


Data=[Data;d];
   
save collect Data

return

d=[str2num(filename(5:8)) MeanDTrpc/1000];   
Data=[Data;d];
   
save collect Data
return
