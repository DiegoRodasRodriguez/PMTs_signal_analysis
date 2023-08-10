% Generic convolution tool for any two distributions in any conjugate space 
% with arbitrary ranges and bins (Some assumptions are made on the response
% function, hence if it represents a generic frequency-domain filter, function
% MTLFilter is better suited)
%
% (Diego Gonzalez-Diaz 19-06-2011)
% Modified (DGD) 10-01-2012
% rating 5*
% 
% Guidelines:
% 1) Function 1 is treated as main function and function 2 is treated as
% a response function. This might have no influence in the result but it is
% strongly advised to follow this convention. This means in particular that
% function 2 must start at 0. ***If an offset is needed for that, it should be
% subtracted at the end.***
% 2) It is assumed that function 2 is properly normalized (i.e., it is a pdf)
% Otherwise MTLconvol will detect it an cure it (and verbose it). If not a pdf
% the result must be multiplied by the function's maximum after convolution.
% 3) For cases where time and frequency domains are combined make sure that
% units are conjugate [t] == 1/[f] (*do not use angular frequency*).
% 4) Bins have to be uniformly spaced, although they might differ in general
% for function 1 and 2.
% 5) Ranges do not need to coincide, but chose reasonable ones in your own benefit.
%
% NOTE: there are some parts that are done with non-optimal function 
% interp1 for simpler code writing. Might be worth improving it for 
% very large arrays.
%
% function [I, Time] = MTLconvol(U1, X1, U2, X2, Mode, Type);
%
% INPUT VARIABLES:
% U1;    Function 1 (complex in general)
% X1;    X-variable 1
% U2;    Normalized Distribution Function 2 (complex in general)
% X2;    X-variable 2
% Mode;  The input variables can be time-freq, freq-freq, freq-time,
%        time-time. It is specified by this variable as:
%        'tf', 'ff', 'ft', 'tt'.
% Type=='FFT', 'Polynomial'. Default 'FFT'.
% verbose; level of verbosity. Default 0.
%
% OUTPUT VARIABLES:
% I;     Convoluted function
% Time;  Time
%
% AUXILIARY FUNCTIONS:
%
%function [F_sym, X_sym, FWasSym] = MTLsymmetrizeFFT(F, X)
%function [F_asym, X_asym]        = MTLasymmetrizeFFT(F, X)
%
% AUXILIARY SUBFUNCTIONS:
%
%function [U1_, X1_, U2_, X2_, RangesAreTheSame] = ...
%   MTLmatchArrays(U1, X1, U2, X2)
%

function [I, Time] = MTLconvol(U1, X1, U2, X2, Mode, Type, verbose)

if nargin<7,  verbose = 1; end;
if nargin<6,  Type='FFT';  end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% First run some consistency checks to increase robustness %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Convert all arrays to rows for simplified treatment
sizeU1 = size(U1); if(sizeU1(1)~=1), U1 = U1'; end
sizeX1 = size(X1); if(sizeX1(1)~=1), X1 = X1'; end
sizeU2 = size(U2); if(sizeU2(1)~=1), U2 = U2'; end
sizeX2 = size(X2); if(sizeX2(1)~=1), X2 = X2'; end

%End with trivial output if any of the two input arrays is zero
if(isempty(U1(U1~=0)) || isempty(U2(U2~=0)))
    I    = zeros(size(U1));
    Time = X1;
    return;
end

%Check that X-arrays have uniform bin spacing
isIsoX1 = MTLcheckThatArrayIsIsoSpaced(X1);
isIsoX2 = MTLcheckThatArrayIsIsoSpaced(X2);

if(~isIsoX1 || ~isIsoX2), error('***X1, X2 arrays not iso-spaced!***'); end

%Convert first to time domain, symmetrizing the response functions if needed
if Mode(1)=='f'
    if(X1(1)~=0) 
        warning('*** DC value not specified, include unless it comes from data ***'); %#ok<WNTAG>
    end  
    [U1_ X1_] = MTLsymmetrizeFFT(U1,X1);
    U1 = ifft(U1_);
    X1 = (0:(length(X1_)-1))/max(X1_);
end
if Mode(2)=='f'
    if(X2(1)~=0) 
        warning('*** DC value not specified, include unless it comes from data ***'); %#ok<WNTAG>
    end  
    [U2_ X2_] = MTLsymmetrizeFFT(U2,X2);
    U2 = ifft(U2_);
    X2 = (0:(length(X2_)-1))/max(X2_);
end

%Match ranges
[U1_, X1_, U2_, ~, RangesAreTheSame] = MTLmatchArrays(U1, X1, U2, X2);
stept = X1_(2) - X1_(1);

%Re-define h-response in order to be normalized for numerical treatment:
U2_  = stept*U2_;

%Check that the discrete sum of the response function is one

if((abs(sum(U2_)-1) > 0.1))
    U2_ = U2_/sum(abs(U2_));
    if(verbose)
        warning('Integral of response function not 1, within more than 10%'); %#ok<WNTAG>
        disp('Automatic re-normalization performed');
    end
end
%%%%%%%%%%%%%%%%%%%%%%
%%% Do convolution %%%
%%%%%%%%%%%%%%%%%%%%%%

if(strcmpi(Type,'FFT')==1)
    I        = real(ifft(fft(U2_).*fft(U1_)));
    Time     = X1_;
elseif(strcmpi(Type,'Polynomial')==1)
    I        = conv(U2_, U1_);
    I        = I(1:(length(X1_)));
    Time     = X1_;
else
    error('***Incorrect convolution type in function MTLconvol***');
end

%Re-calculate time range set to fit to the one of the main function
if(~RangesAreTheSame)
    I    = interp1(X1_, I, X1);
    I    = I(~isnan(I));
    Time = X1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Sub-function match ranges of different function arrays %%%%%%

function [U1_, X1_, U2_, X2_, RangesAreTheSame] = ...
    MTLmatchArrays(U1, X1, U2, X2)

U1_=U1; X1_=X1; U2_=U2; X2_=X2;
RangesAreTheSame = 0;

if(isequal(X1,X2))  %X-arrays are identical (skip the rest)
    RangesAreTheSame = 1; 
    return; 
else                %X-arrays are not identical (match)
    
    stept = X1(2)-X1(1);
    
    if(max(X1)>max(X2)), Xmax = max(X1);
    else                 Xmax = max(X2);
    end
    if(min(X1)<min(X2)), Xmin = min(X1);
    else                 Xmin = min(X2);
    end
    %Numerical correction to keep spacing
    Nsteps = ceil((Xmax-Xmin)/stept);
    Xmax   = Xmin + Nsteps*stept;
    Xrange = Xmin:stept:Xmax;
    
    U2_ = interp1(X2, U2, Xrange);
    U2_(isnan(U2_)) = 0;
    
    if(Xmax ~= max(X1) || Xmin ~= min(X1)) %Might be needed if change in Xmax
        U1_ = interp1(X1, U1, Xrange);  %Not the fastest algorithm for this
        U1_(isnan(U1_)) = 0;
    end
    
    X1_ = Xrange;
    X2_ = Xrange;
    
end
