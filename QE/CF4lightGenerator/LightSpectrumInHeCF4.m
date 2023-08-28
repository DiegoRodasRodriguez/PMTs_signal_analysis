%This is a pretty old He/CF4 generator that needs of some renovation...
%in particular including the references is clearly important.
%Use Vienna work for reference!.
%Obtain interpolated spectrum for HeCF4 secondary
%scintillation in GEMs.
%The only free variable is the fraction of CF4,
%but it could be possibly done more general to include
%pressure-dependence, field, S1 or S2.
%Also wavelength range is important
%FIXME: detail the reference for each spectrum
%At the moment it only works for discrete values of Lambda_...

function [Y, X] = LightSpectrumInHeCF4(fCF4, Lambda_)

YieldInterp   = zeros(size(Lambda_));

load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_PureCF4_1bar_GEM_NORM;
LambdaPure = Lambda; YieldPure = Yield;
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_HeCF4_60_40_1bar_GEM_NORM;
LambdaHe60 = Lambda; YieldHe60 = Yield; 
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_PureHe_1bar_GEM_NORM;
LambdaPureHe = Lambda; YieldPureHe = Yield;

for k=1:length(Lambda_)
    YieldPure_   = YieldPure  (LambdaPure   == Lambda_(k)); if(isempty(YieldPure_)),   YieldPure_= -1;   end
    YieldHe60_   = YieldHe60  (LambdaHe60   == Lambda_(k)); if(isempty(YieldHe60_)),   YieldHe60_= -1;   end
    YieldPureHe_ = YieldPureHe(LambdaPureHe == Lambda_(k)); if(isempty(YieldPureHe_)), YieldPureHe_= -1; end
    
    YieldBin   = [YieldPure_, YieldHe60_, YieldPureHe_];
    ConcBin    = [1, 0.4, 0];
    YieldBin_  = YieldBin(YieldBin>=0);
    ConcBin_   = ConcBin (YieldBin>=0);
    if(isempty(YieldBin_)),       YieldInterp(k)=0;
    elseif(length(YieldBin_)==1), YieldInterp(k)=YieldBin_;
    else YieldInterp(k) = interp1(ConcBin_, YieldBin_, fCF4, 'linear', 'extrap');
    end
end
YieldInterp(YieldInterp<0) = 0;

Y=YieldInterp;
X=Lambda_;

end