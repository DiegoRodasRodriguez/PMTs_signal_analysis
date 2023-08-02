function [Mean, Sigma, Max]=stdfitGraph_D(X,Y, Nsigmas, sigmaY, titulo)

global Xdata Ydata sigmaData;

plotit=0;
if nargin==5,  plotit=1;             end
if nargin==2,  Nsigmas=1.5;          end
if nargin<=3,  sigmaY=ones(size(X)); end

ULimit  = max(X); 
LLimit  = min(X);
AbsY    = abs(Y);

Xdata = X; Ydata = Y; sigmaData = sigmaY;

par = [mean(X), ULimit - LLimit, Y(AbsY==max(AbsY))];

lb  = [LLimit, 0,               min(0,2*Y(AbsY==max(AbsY)))];
ub  = [ULimit, ULimit - LLimit, max(0,2*Y(AbsY==max(AbsY)))];
    
for i=1:3    
    
    options = optimset('TolX',1e-9,'TolFun',1e-13,'MaxFunEvals', 1000*length(par),'MaxIter',2000);
    
    par     = lsqnonlin(@gChi2,par,lb,ub,options);
    
    Mean    = par(1);
    Sigma   = par(2);
    Max     = par(3);
    Xmax    = min(max(X), Mean + Sigma*Nsigmas);
    Xmin    = max(min(X), Mean - Sigma*Nsigmas);
    
    Xdata     = X(X>=Xmin & X<=Xmax);
    Ydata     = Y(X>=Xmin & X<=Xmax);
    sigmaData = sigmaY(X>=Xmin & X<=Xmax);
    
end

if plotit 
    figure; plot(Xdata, Ydata, 'b'); hold on;
    plot(Xdata, par(3)*exp(-(Xdata-par(1)).^2/(2*par(2)^2)), 'g');
	yaxis(min(Y), max(Y));
	xaxis(min(X), max(X));
	title(titulo);
	%drawnow;
end

return