function [stdev, chi2, par]=sigma(x,n)
% returns a gauss fit to the histogram x ,n
% par(2)=sigma; par(1)=mean; par(3)=const

idx=find(n>1);
x=x(idx);
y=log(n(idx));
x=x(:);
y=y(:);

[p,S] = polyfit(x,y,2);

p

A=[x.^2 x x*0+1];
p=lscov(A,y,diag(y))

max(A*p-y)
%keyboard

stdev=sqrt(-1/p(1)/2);
mean=-p(2)/2/p(1);
const=exp((-p(2)^2+4*p(3)*p(1))/p(1)/4);

par=[mean, stdev,const];

ndf=length(n)-3;

chi2 = sum((exp(y) - g(x,par)).^2 ./ exp(y))/ndf;


