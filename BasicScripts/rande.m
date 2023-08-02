function [out]=rande(n,a)

% Generate a negative exp distribution y=exp(-ax)
% rande(number of events,coef of exp)

out=real(log(-rand(n,1)/a)*(-1/a));

return;
