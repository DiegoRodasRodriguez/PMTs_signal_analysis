function plotbar(x,y,Sigmas)

%   Dibuja x frente a y con las barrar de eror 
%   de sigma
%
%


figure
plot(x,y,'xb');lastline('markersize',6);

hold on;
plot(x,y,'-b');
plot((x(:)*[1 1])',(y(:)*[1 1]+Sigmas(:)*[-1 1])','-r');
return