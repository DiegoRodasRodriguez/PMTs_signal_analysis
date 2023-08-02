function [parout,chisq]=chisqmin(par, xin, yin, sigmain)

global x;
global y;
global sigma;
global ndf;

clear parout;

x = xin;
y = yin;
if nargin == 4
  sigma = sigmain;
else
  sigma = sqrt(y);
end;

ndf=length(find (sigma > 0));
%parout=fmins('gx',par);
parout=fminsearch('gx',par);
chisq = gx(parout);
