%
% Simple function to cross-check that array is linearly iso-spaced
% DGD (29-01-2013)
%
% function isIso = MTLcheckThatArrayIsIsoSpaced(X)
%
% INPUT VARIABLES:
% X;     array to be checked
%
% OUTPUT VARIABLES:
% isIso; is 1, it is iso-spaced
%

function isIso = MTLcheckThatArrayIsIsoSpaced(X)
isIso = 0;
stepX = min(diff(X));
if((abs((min(diff(X))-max(diff(X))))/stepX < 1e-10)), isIso =1; end
end