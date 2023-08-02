function [Mean, Sigma, Max, Tails, Tails300]=stdfit_new_D(Tcorr, n_sigmas, titulo, low_T, up_T, step_size)

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
%NOTE: The fitting function hfitg_D is a disgrace
nbins=ceil((up_T-low_T)/step_size);
par =hfitg_D(Tcorr,nbins,0,low_T,up_T);
Mean=par(1);Sigma=par(2);Max=par(3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Exception to avoid breakdown of iterative fit
if(Mean<low_T || Mean>up_T),                                           Mean  = (low_T + up_T)/2;      end
if(Sigma>(up_T-low_T)/sqrt(12) || Sigma<0),                            Sigma = (up_T-low_T)/sqrt(12); end
if(Mean<low_T || Mean>up_T || Sigma>(up_T-low_T)/sqrt(12) || Sigma<0), Max   = Max_T;                 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Limit=Sigma*n_sigmas;    
% Fit in +-n_sigmas around the maximum
eps=0.05;  % Convergency criteria
nmax=30;

for i=1:nmax
%     low_T = low_T
%     up_T  = up_T
%     Mean  = Mean
%     Sigma = Sigma
%     Max   = Max
%     i=i
%     Limit = Limit
    Mean_last=Mean;
    [N,X]=hist1D(Tcorr,-Limit+Mean:step_size:Limit+Mean);
    par =hfitg_D(X,N,0);
    Mean=par(1) ; Sigma=par(2);Max=par(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Exception to avoid breakdown of iterative fit
    if(Mean<low_T || Mean>up_T),                                           Mean  = (low_T + up_T)/2;      end
    if(Sigma>(up_T-low_T)/sqrt(12) || Sigma<0),                            Sigma = (up_T-low_T)/sqrt(12); end
    if(Mean<low_T || Mean>up_T || Sigma>(up_T-low_T)/sqrt(12) || Sigma<0), Max   = Max_T;                 end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Limit=Sigma*n_sigmas;
    ratio=abs((Mean_last-Mean)*2/(Mean_last + Mean));
    if ratio<eps, break;end;    
end

Tcorr=Tcorr-Mean;
TailsRight=length(find(Tcorr>3*Sigma))/length(Tcorr);
TailsLeft=length(find(Tcorr<-3*Sigma))/length(Tcorr);
Tails=TailsLeft+TailsRight;
Tails300=length(find(abs(Tcorr)>300))/length(Tcorr);

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
    xlabel('Time difference [ps]','FontSize',FontSize+2)
    ylabel(['Events/' num2str(step_size) ' ps'],'FontSize',FontSize+2)
    
    text(min(xaxis)+diff(xaxis)/20,max(yaxis)/2,['Sigma=' num2str(Sigma)],'FontSize',FontSize);
    text(min(xaxis)+diff(xaxis)/20,max(yaxis)/3,['3-Sigma Tails=' num2str(Tails*100) '%'],'FontSize',FontSize);
    text(min(xaxis)+diff(xaxis)/20,max(yaxis)/5,['300ps Tails=' num2str(Tails300*100) '%'],'FontSize',FontSize);
    text(min(xaxis)+diff(xaxis)/20,max(yaxis)/10,['Events=' num2str(length(Tcorr))],'FontSize',FontSize);
    title(titulo);
    drawnow
    
end

Mean=Mean+Max_T+Mean_T;

return
