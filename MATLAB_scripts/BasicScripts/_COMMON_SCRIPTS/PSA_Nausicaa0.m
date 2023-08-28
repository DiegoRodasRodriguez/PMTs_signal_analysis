function [Int, Amp, tMin, tMax, FWHM, meanT, trise, tfall, rmsOUT, trise20, tfall20, Ped, assy, stdT] = PSA_Nausicaa0(F, time, thres, trigTime, tit)

%Calculate charge and amplitude
tMin  = 0.85*max(time(F<thres & time<trigTime)); if(isempty(tMin)), tMin=time(1);            end
tMax  = 1.3*min(time(F<thres & time>trigTime));  if(isempty(tMax)), tMax=time(length(time)); end
Int   = sum(F(time>=tMin & time<=tMax));
Amp   = max(F(time>=tMin & time<=tMax));

%Calculate FWHM and meanTime
tAtmax = mean(time(F==Amp));  
FWHML  = max(time(F<Amp/2 & time<tAtmax)); if(isempty(FWHML)), FWHML=0; end
FWHMR  = min(time(F<Amp/2 & time>tAtmax)); if(isempty(FWHMR)), FWHMR=0; end

FWHM   = FWHMR - FWHML;
%meanT  = sum(F(time>=tMin & time<=tMax).*time(time>=tMin & time<=tMax))/sum(F(time>=tMin & time<=tMax));
%meanT2 = sum(F(time>=tMin & time<=tMax).*time(time>=tMin & time<=tMax).*time(time>=tMin & time<=tMax))/sum(F(time>=tMin & time<=tMax));
meanT  = sum(abs(F(time>=tMin & time<=tMax)).*time(time>=tMin & time<=tMax))/sum(abs(F(time>=tMin & time<=tMax)));
meanT2 = sum(abs(F(time>=tMin & time<=tMax)).*time(time>=tMin & time<=tMax).*time(time>=tMin & time<=tMax))/sum(abs(F(time>=tMin & time<=tMax)));
stdT   = sqrt(meanT2 - meanT.^2);
stdT   = real(stdT);

%Calculate rise and fall time (NOTE: not from 0.1 to 0.9 as usual, but peak to noise)
trise  = mean(time(F==Amp)) - tMin;
tfall  = tMax - mean(time(F==Amp));

%Calculate rise and fall time (NOTE: not from 0.1 to 0.9 as usual), but peak to 0.2
rise20left  = max(time(F<Amp*0.2 & time<tAtmax)); if(isempty(rise20left)),  rise20left=0;  end
rise20right = min(time(F<Amp*0.2 & time>tAtmax)); if(isempty(rise20right)), rise20right=0; end

trise20 = tAtmax - rise20left;
tfall20 = rise20right - tAtmax;

%RMS outside window (can be used for removing oscillating events, very noisy, or simply strange events)
rmsOUT = std(F(time<tMin | time>tMax));

AmpInsideWindow  = F(time>=tMin & time<=tMax);
AmpOutsideWindow = F(time<tMin | time>tMax);

if(length(AmpOutsideWindow)>=length(AmpInsideWindow)), Ped = sum(AmpOutsideWindow(1:length(AmpInsideWindow)));
else, Ped = -999;
end

assy = -(FWHMR - meanT)/(FWHMR - FWHML);

if nargin>4
    figure; hold on;
    plot(time,F,'-');
    line([trigTime, trigTime], [min(F),max(F)]); lastline('color','r'); lastline('LineStyle','--');
    line([tMin, tMin], [min(F),max(F)]); lastline('color','r');
    line([tMax, tMax], [min(F),max(F)]); lastline('color','r');    
    line([tAtmax, tAtmax], [min(F),max(F)]); lastline('color','c'); lastline('LineStyle','--');
    line([FWHML, FWHMR], [max(F)/2,max(F)/2]);      lastline('color','k');    
    line([rise20left, rise20right], [max(F)/4,max(F)/4]);      lastline('color','k');
    xlabel('time [us]'); ylabel('A[mV]'); title(tit);
    box on;
    disp(' ');
    disp('************ PULSE parameters **********');
    disp(['integral= ', num2str(Int), ' Amp= ', num2str(Amp), '[mV], tmin= ', num2str(tMin), ' us, tmax= ', num2str(tMax), ' us ']);
    disp(' ');
    disp(['FWHM= ', num2str(FWHM), ' us, stdT=' num2str(stdT), ' us, trise= ', num2str(trise), ' us, tfall= ', num2str(tfall), ' us ']);
    disp(' ');    
    disp(['assy= ', num2str(assy)]);
    disp(' ');
        
end
