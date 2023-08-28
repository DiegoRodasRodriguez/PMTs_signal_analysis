function [parout,chisq]=hfitg_D(x,n,~,miny,maxy)
%Gaussian Fit

if all(size(x)==size(n))
	HistIsGiven=1;
else
	HistIsGiven=0;
end

if (nargin == 5)  % trim x
  x=x((x>=miny & x<=maxy));
else
	miny=min(min(x));
	maxy=max(max(x));
end

if HistIsGiven,	nx=x; ny=n;
else            [ny,nx]=hist(x,n,miny,maxy); end

par(1)=mean(x);
par(2)=std(x);
par(3)=max(ny);
[parout,chisq]=chisqmin(par,nx,ny);
