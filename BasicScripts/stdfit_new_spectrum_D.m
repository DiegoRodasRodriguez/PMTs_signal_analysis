%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function for recursive fit and peak finding
%DGD (30/11/2015)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Mean, Sigma, Max, Nevts]=stdfit_new_spectrum_D(Tcorr, n_sigmas, titulo, low_T, up_T, step_size)

FontSize=8;
plotit=0;

if nargin>2, plotit=1;end

Mean_T = mean(Tcorr);
Tcorr  = Tcorr-Mean_T;

if nargin<=3
   low_T     = min(Tcorr);
   up_T      = max(Tcorr);
   step_size = (up_T-low_T)/100;
end

[N,X]  = hist1D(Tcorr, low_T:step_size:up_T); 
Max_T  = max(X(N==max(N)));
Tcorr  = Tcorr-Max_T;

%First fit for getting the mean and the sigma
%nbins=ceil((up_T-low_T)/step_size); 
%par =hfitg_D(Tcorr,nbins,0,low_T,up_T);
%Mean=par(1);Sigma=par(2);Max=par(3);
Mean  = mean(Tcorr); Sigma=std(Tcorr);
Limit = Sigma*5;
up_T  = Mean + Limit;
low_T = Mean - Limit;
%nbins=ceil((up_T-low_T)/step_size); 
 
% Fit in +-n_sigmas around the maximum
eps=0.05;  % Convergency criteria
nmax=30;

%iCount = 0;
for i=1:nmax
%    i=i
    Mean_last = Mean;    
    step_size = 2*Limit/100;
    [N,X]=hist1D(Tcorr,-Limit+Mean:step_size:Limit+Mean);
    par =hfitg_D(X,N);
    if(par(2)<0) %Bad Fit (patch to continue)
        Sigma = std(Tcorr(Tcorr>(-Limit+Mean) & Tcorr<(+Limit+Mean)));
        Mean  = mean(Tcorr(Tcorr>(-Limit+Mean) & Tcorr<(+Limit+Mean)));
        Max   = max(N);
        Limit = Sigma*n_sigmas;
        continue;
    else
        Mean=par(1); Sigma=par(2); Max=par(3);
%        disp('in fit');
    end        
    Limit=Sigma*n_sigmas;
    ratio=abs((Mean_last-Mean)*2/(Mean_last + Mean));
    if ratio<eps, break;end;
end

Nevts = Max*(Sigma*sqrt(2*pi))/step_size;

Tcorr      = Tcorr-Mean;

if plotit
    
x=low_T:step_size:up_T-step_size;
[N,X]=hist1D(Tcorr,x);
figure;stairs(X,N); hold on;lastline('color','k')
x2=-n_sigmas*Sigma:5:n_sigmas*Sigma;
plot(x2,g(x2,0,Sigma)*Max,'k')
lastline('linewidth',2);logy;
plot(x,g(x,0,Sigma)*Max,'--k','linewidth',1)
yaxis(1,2*max(N));
set(gca,'FontSize',FontSize)
title(titulo);
xlabel('charge units','FontSize',FontSize+2)
ylabel(['Events/' num2str(step_size) ' ps'],'FontSize',FontSize+2)

text(min(xaxis)+diff(xaxis)/20,max(yaxis)/2,['Sigma=' num2str(Sigma)],'FontSize',FontSize);
text(min(xaxis)+diff(xaxis)/20,max(yaxis)/10,['Events=' num2str(length(Tcorr))],'FontSize',FontSize);
title(titulo);
drawnow

end

Mean=Mean+Max_T+Mean_T;

return
