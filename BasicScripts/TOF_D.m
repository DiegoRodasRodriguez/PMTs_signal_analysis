T_=T; A_=A;
MeanT=mean(T_);RawSigmaT=std(T_);

TimeRange=MeanT+20000*[-1 1];
OldRawSigmaT=1e10;
k=1;
while ((OldRawSigmaT-RawSigmaT)/OldRawSigmaT>.2 | RawSigmaT>1000) & k<30
    k=k+1;
    TimeRange=MeanT+3.0*RawSigmaT*[-1 1];
    idx=find(T_>TimeRange(1) & T_<TimeRange(2));
    T_=T_(idx); A_=A_(idx);
    MeanT=mean(T_);
    OldRawSigmaT=RawSigmaT;
    RawSigmaT=std(T_);
end

TimeRange=MeanT+3*RawSigmaT*[-2.0 4];

%TimeRange=MeanT+1500*[-1 1];%%

idx=find(T>TimeRange(1) & T<TimeRange(2));
T_=T(idx); A_=A(idx);

RawSigmaT=std(T_);


% Mean values ---------------------------------------
MeanT=mean(T_);


%MeanT=max(X(find(N==max(N)))); %%
%MeanT = MeanT

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
[As_,I]=sort(A_);
Ts_=T_(I);
clear I
lAs_=length(As_);
coefs=[];
% subdivide first bin
Abnd=floor([linspace(1,lAs_/NrBounds*3/4,3) linspace(lAs_/NrBounds,lAs_,NrBounds-3)]); 

for i=2:NrBounds
    Tt_=Ts_(Abnd(i-1)+1:Abnd(i));
    Aa_=As_(Abnd(i-1)+1:Abnd(i));
    MeansT_(i-1)=mean(Tt_);
    Centers(i-1)=mean(Aa_);  % centers of the A-bins
    
    p=[0 polyfit(Aa_,Tt_,1)];
    Tt_=(Tt_-polyval(p,Aa_));
    
    coefs=[coefs; p];
    
    Sigmas=zeros(1,NrBins)+1;
    MeansT__=zeros(1,NrBins);
    
    A__=[A__; Aa_]; 
    T__=[T__; Tt_]; % corrected T_
    
end

%Sigmas=sqrt(Sigmas.^2-SigmaStart^2);

p=polyfit(Centers(2:6),MeansT_(2:6),2); % special parabola for extrapolation at small A
coefs(1,:)=p;   % insert parabola
TApp=mkpp(As_(Abnd(:))',coefs(:)); %standard fit

Tcorr=T-myppval(TApp,A);
RawMeanT=mean(MeansT_);


% T-A plots ------------------------------------------------------------------

figure 

NrScat=floor(rand(1,500)*(length(T__)-1)+1);
subplot(2,2,1);
plot(A_(NrScat),T_(NrScat),'.r'); hold on; 
plot(Centers,MeansT_,'ob');
axis([0 1000 TimeRange]);
plot(As_(Abnd),myppval(TApp,As_(Abnd)),'+b');
drawnow

subplot(2,2,2); 
plot(A__(NrScat),T__(NrScat),'.'); hold on;
plot(Centers,MeansT__,'or');
plot((Centers(:)*[1 1])',(MeansT__(:)*[1 1]+Sigmas(:)*[-1 1])','-r');
axis([0 2000 -1000 1000]);
drawnow

subplot(2,2,3);
[N,X]=histf(T_-MeanT,linspace(-1500,1500,150)); stairs(X,N+1);logy; 
xaxis(-1500,1500);
drawnow

subplot(2,2,4);
[N,X]=histf(T__,linspace(-1500,1500,150)); stairs(X,N+1);logy;
xaxis(-1500,1500);
drawnow
