%This is a pretty old Ar/CF4 generator that needs of some renovation...
%in particular including the references is clearly important.
%Use Vienna work for reference!.
%Obtain interpolated spectrum for ArCF4 secondary
%scintillation in GEMs.
%The only free variable is the fraction of CF4,
%but it could be possibly done more general to include
%pressure-dependence, field, S1 or S2.
%Also wavelength range is important
%FIXME: detail the reference for each spectrum
%At the moment it only works for discrete values of Lambda_ at 1nm
%NOTE: 16/01/2023. Find implementation a bit weird :)... the more or
%less interesting key is that it gives an answer for every concentration

function [Y, X] = LightSpectrumInArCF4_sec(fCF4, Lambda_)

YieldInterp   = zeros(size(Lambda_));

load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_PureCF4_1bar_GEM_NORM;
LambdaPure = Lambda; YieldPure = Yield;
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_ArCF4_95_5_1bar_GEM_NORM;
LambdaAr95 = Lambda; YieldAr95 = Yield;
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_ArCF4_90_10_1bar_GEM_UVcut_NORM;
LambdaAr90 = Lambda; YieldAr90 = Yield;
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_ArCF4_33_67_1bar_GEM_NORM;
LambdaAr33 = Lambda; YieldAr33 = Yield;
load E:\HOME_RareEventsGroup\Diego\QE_Calc\CF4lightGenerator\S2_PureAr_1bar_GEM_NORM;
LambdaPureAr = Lambda; YieldPureAr = Yield;

for k=1:length(Lambda_)
    YieldPure_   = YieldPure  (LambdaPure   == Lambda_(k)); if(isempty(YieldPure_)),   YieldPure_= -1;   end
    YieldAr95_   = YieldAr95  (LambdaAr95   == Lambda_(k)); if(isempty(YieldAr95_)),   YieldAr95_= -1;   end
    YieldAr90_   = YieldAr90  (LambdaAr90   == Lambda_(k)); if(isempty(YieldAr90_)),   YieldAr90_= -1;   end
    YieldAr33_   = YieldAr33  (LambdaAr33   == Lambda_(k)); if(isempty(YieldAr33_)),   YieldAr33_= -1;   end
    YieldPureAr_ = YieldPureAr(LambdaPureAr == Lambda_(k)); if(isempty(YieldPureAr_)), YieldPureAr_= -1; end
    
    YieldBin   = [YieldPure_, YieldAr33_, YieldAr90_, YieldAr95_, YieldPureAr_];
    ConcBin    = [1, 0.67, 0.1, 0.05, 0];
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