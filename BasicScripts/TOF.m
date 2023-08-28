idx=find(T<(2048) & T>0); % Discard overflow channels
T=T*TDC2binSize;
T_=T(idx); A_=A(idx);
%idx=find(A_>200); % Discard amplitude underflows
%T_=T_(idx); A_=A_(idx);

% reject far away events--------------------------------------------------
% forbiden in case of tails

MeanT=mean(T_);RawSigmaT=std(T_)
TimeRange=MeanT+20000*[-1 1];
OldRawSigmaT=1e10;
k=1;
%while ((OldRawSigmaT-RawSigmaT)/OldRawSigmaT>.2 | RawSigmaT>1000) & k<30
while ((OldRawSigmaT-RawSigmaT)/OldRawSigmaT>.2 | RawSigmaT>1000) & k<30
    k=k+1;
    %	TimeRange=MeanT+4*RawSigmaT*[-1 1];
    TimeRange=MeanT+2.5*RawSigmaT*[-1 1];
    idx=find(T_>TimeRange(1) & T_<TimeRange(2));
    T_=T_(idx); A_=A_(idx);
    MeanT=mean(T_);
    OldRawSigmaT=RawSigmaT;
    RawSigmaT=std(T_);
end

TimeRange=MeanT+3*RawSigmaT*[-1 4];
%TimeRange=MeanT+3*RawSigmaT*[-2 3]; % standard
%TimeRange=MeanT+3*RawSigmaT*[-2 6]; % mostly for electronic calib with varamp pulser
idx=find(T>TimeRange(1) & T<TimeRange(2));
T_=T(idx); A_=A(idx);

RawSigmaT=std(T_)
% all T_s are now inside time range
% clear T A

%------ cut lower amplitudes ---------------------------------------

if 0
    
    if all (A_<10)
        [N,X]=histf(A_,0:.05:3); 
        Llim=X(min(find(N/max(N)>0.05)));
        figure; stairs(X,N/max(N)); logy
        (length(find(A_>Llim)))/length(A_)
        [A_,T_]=cut(find(A_>Llim),A_,T_);
        [N,X]=histf(A_,0:.05:3); hold on; stairs(X,N/max(N)); lastline('color','b')
    else
        [N,X]=histf(A_,0:10:2048); 
        Llim=X(min(find(N/max(N)>0.1)));
        figure; stairs(X,N/max(N)); logy
        (length(find(A_>Llim)))/length(A_)
        [A_,T_]=cut(find(A_>Llim),A_,T_);
        [N,X]=histf(A_,0:10:2048); hold on; stairs(X,N/max(N)); lastline('color','b')
    end
    
end

% Mean values ---------------------------------------

MeanT=mean(T_);
if all(A_<10)
    I=find(A_<=2); pMean=polyfit(A_(I),T_(I),1);
else
    pMean=0;
end

% slice the amplitude spectrum and correct ---------------------------------------
NrBounds=13;	% 23 for many events
%NrBounds=4;	

NrBins=NrBounds-1;
Sigmas_=zeros(1,NrBins);
Centers=zeros(1,NrBins);

A__=[];T__=[];
%ACorr=[];TCorr=[];

[As_,I]=sort(A_);
Ts_=T_(I);
clear I
lAs_=length(As_);
coefs=[];

%Eff_TOF_Software=length(Ts_)/Events*100

%Abnd=[1 min(find(As_>1.2)) linspace(min(find(As_>1.4)),lAs_,NrBounds-2)];
%Abnd=[linspace(1,lAs_/NrBounds*10/11,10) linspace(lAs_/NrBounds,lAs_,NrBounds-10)]; % subdivide first bin
Abnd=[linspace(1,lAs_/NrBounds*3/4,3) linspace(lAs_/NrBounds,lAs_,NrBounds-3)]; % subdivide first bin
%Abnd=floor([linspace(1,lAs_,NrBounds)]); 

for i=2:NrBounds
    Tt_=Ts_(Abnd(i-1)+1:Abnd(i));
    Aa_=As_(Abnd(i-1)+1:Abnd(i));
    MeansT_(i-1)=mean(Tt_);
    Centers(i-1)=mean(Aa_);  % centers of the A-bins
    if all(A_<10)   % log scale
        % 		p=polyfit(Aa_,Tt_,2);
        p=[0 polyfit(Aa_,Tt_,1)];
    else
        %		p=[0  0 polyfit(Aa_,Tt_,2)];
        p=[0 polyfit(Aa_,Tt_,1)];
    end
    Tt_=(Tt_-polyval(p,Aa_));
    
    coefs=[coefs; p];
    if 0 % sigma measurement in each bin (0=skip)
        %	[Mean, Sigma, Max]=stdfit(Tt_,'test');
        [Mean, Sigma, Max]=stdfit(Tt_);
        Sigmas(i)=Sigma;
        MeansT__(i)=Mean;
    else
        Sigmas=zeros(1,NrBins)+1;
        MeansT__=zeros(1,NrBins);
    end
    A__=[A__; Aa_]; 
    T__=[T__; Tt_]; % corrected T_
    
end

Sigmas=sqrt(Sigmas.^2-SigmaStart^2);

p=polyfit(Centers(2:6),MeansT_(2:6),2); % special parabola for extrapolation at small A
coefs(1,:)=p;   % insert parabola
TApp=mkpp(As_(Abnd(:))',coefs(:)); %standard fit

Tcorr=T-myppval(TApp,A);

findobj('UserData','rmsTOF');
if length(ans)==0
    figure
    marker='o';
else
    axes(ans);
    marker='*';
    text(800,500,'RPC - o');
    text(800,450,'PPC - *');
end

plot(Centers,Sigmas,marker); hold on
if all(A<10)
    axis([1 log10(2400) 0 600]); grid on;
else
    axis([0 2400 0 600]); grid on;
end
xlabel('Amplitude');ylabel('sigma TOF');
set(gca,'UserData','rmsTOF');
drawnow

CorrRmsT=sqrt(std(T__).^2-SigmaStart^2);

RawMeanT=mean(MeansT_)


if 0
    %------------- new fit
    
    x0=[70 560*50 -3*50]
    %I=find(T_>MeanT);
    I=find(A_<80);
    
    T_=T_(I);A_=A_(I);
    x=fmins('ATchi2',x0,foptions,[],A_,T_)
    A__=A_;T__=T_-TAfun(x,A_);
    %I=find(T>MeanT);
    I=find(A<80);
    Tcorr=T(I)-TAfun(x,A(I));
end

% T-A plots ------------------------------------------------------------------
figure 

NrScat=rand(1,500)*(length(T__)-1)+1;
%NrScat=1:length(T__);
subplot(2,2,1);
%plot(A,T,'.c'); hold on
plot(A_(NrScat),T_(NrScat),'.r'); hold on; 
ylabel('Measured time(ps)'); xlabel('log10(Charge/bin)')
plot(Centers,MeansT_,'ob');
%plot((Centers(:)*[1 1])',(MeansT_(:)*[1 1]+Sigmas(:)*[-1 1])','-r');

if length(findstr(titulo,'RPC'))
    if all(A<10)
        %		axis([1 log10(2400) mean(TimeRange)-1000 mean(TimeRange)+1000]);  % for module
        axis([floor(min(A)) log10(2048) TimeRange]);  % for module
    else
        axis([0 1000 TimeRange]);
    end
else
    if all(A<10)
        %		axis([1 log10(2400) mean(TimeRange)-1000 mean(TimeRange)+1000]);
        axis([floor(min(A)) log10(2048) TimeRange]);  % for module
    else
        axis([0 1000 TimeRange]);
    end
    
end

if all(A<10)
    plot(log10(MeanA-min(Arpc)+1),MeanT,'+k'); % mean values of A and T
    plot(linspace(min(A),log10(2048),100),myppval(TApp,linspace(min(A),log10(2048),100)),'b'); % plot the interpolant
    plot([1,2],polyval(pMean,[1,2]),'-g') % mean slope
else
    plot(MeanA,MeanT,'+k'); % mean values of A and T
    a=[min(Arpc):100, 110:10:2400];
    plot(a,myppval(TApp,a),'b'); % plot the interpolant
    %	plot(1:10:2400,TAfun(x0,1:10:2400),'b'); % plot the interpolant
    %	plot(1:10:2400,TAfun(x,1:10:2400),'g'); % plot the interpolant
end
plot(As_(Abnd),myppval(TApp,As_(Abnd)),'+b');

title(titulo);
drawnow

subplot(2,2,2); 
plot(A__(NrScat),T__(NrScat),'.'); hold on; ylabel('Tempo corrigido (ps)');xlabel('log(Amplitude)')
plot(Centers,MeansT__,'or');
plot((Centers(:)*[1 1])',(MeansT__(:)*[1 1]+Sigmas(:)*[-1 1])','-r');
if length(findstr(titulo,'RPC'))
    if all(A<10)
        axis([1 log10(2400) -1000 1000]);  % for module
    else
        axis([0 2400 -1000 1000]);
    end
else
    if all(A<10)
        axis([1 log10(2400) -2000 2000]);
    else
        axis([0 2400 -2000 2000]);
    end
end
drawnow

subplot(2,2,3);
%[N,X]=histf(T_-MeanT,linspace(-CorrRmsT,CorrRmsT,100)*10); stairs(X,N+1);logy; 
[N,X]=histf(T_-MeanT,linspace(-3000,3000,100)); stairs(X,N+1);logy; 
xaxis(-3000,3000);
xlabel('Tempo medido (ps)');
drawnow

subplot(2,2,4);
%[N,X]=histf(T__,linspace(-CorrRmsT,CorrRmsT,100)*10); stairs(X,N+1);logy;
[N,X]=histf(T__,linspace(-3000,3000,100)); stairs(X,N+1);logy;
xaxis(-3000,3000);
xlabel('Tempo corrigido (ps)');
drawnow

% Residual correlations  -------------------------------------------------------

I=find(isnan(Tcorr));Tcorr(I)=0*I;

if 0
    
    disp('---- Residual DTrpc correlation  -----'); % for long counter only
    Tcorr=uncorr(Tcorr,DTrpc);
    
    %disp('---- Residual Tstart correlation  -----');
    %Tcorr=uncorr(Tcorr,Tstart);
    
    %disp('---- Residual Amplitude correlation  -----');
    %Tcorr=uncorr(Tcorr,Arpc);
    
    %disp('---- Residual Baseline correlation  -----');
    %Tcorr=uncorr(Tcorr,Ch22);
    
    %disp('---- Residual scintilator amplitude correlation  -----');
    %Tcorr=uncorr(Tcorr,Aap1);
    
    disp('---- END - Residual DTrpc correlation  -----'); % for long counter only
    
end

% TOF sigma -------------------------------------------------------

[Mean, Sigma, Max]=stdfit2(Tcorr,titulo);Sigma=Sigma
CorrSigmaTOF=sqrt(Sigma^2-SigmaStart^2);

Tcorr=Tcorr-Mean;

end	% Varamp

return


% TOF sigma vs TOF eficiency -------------------------------------------------------

N=length(T__);
Npts=1;
Apts=linspace(1,N/10,Npts); 

TOFeffic=zeros(1,Npts);
TOFrms=zeros(1,Npts);
TOFsigma=zeros(1,Npts);

for i=1:Npts
    Tlocal=T__((Apts(i))+7:N);
    TOFrms(i)=std(Tlocal);
    pars=hfitg(Tlocal,20);	TOFsigma(i)=pars(2);
    TOFeffic(i)=(N-Apts(i))/Events*100;
end

TOFrms=sqrt(TOFrms.^2-SigmaStart^2);
TOFsigma=sqrt(TOFsigma.^2-SigmaStart^2);

axe=findobj('UserData','TOFplot');
if length(ans)==0
    figure
    marker='o';
    set(gca,'UserData','TOFplot');
else
    %	figure(figure(get(axe,'Parent')));
    axes(axe);
    marker='*';
    text(94,140,'RPC - o');
    text(94,120,'PPC - *');
    hold on;
end

plot(TOFeffic,TOFrms,marker); hold on
plot(TOFeffic,TOFsigma,marker); 
title(['TOF rms vs TOF eficiency '])
xlabel('eficiency(%) (from offline amplitude cut)');ylabel ('TOF rms,sigma (ps)');
xaxis(85,100); 
%yaxis(100,300);
yaxis(100,500);

grid on
drawnow
set(gca,'UserData','TOFplot');

