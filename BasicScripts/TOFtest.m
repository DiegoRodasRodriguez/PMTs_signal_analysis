idx=find(T<(2048) & T>30); % Discard overflow channels.

%T=(-T+2*Mean1)*1000/abs(Mean1-Mean2); % if time tail goes to shorter times
%T=T*1000/abs(Mean1-Mean2);% if time tail goes to longer times # use for anarun only

T=T*50;

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
%TimeRange=MeanT+3*RawSigmaT*[-1 2];
TimeRange=MeanT+3*RawSigmaT*[-2 3];
%TimeRange=MeanT+3*RawSigmaT*[-2 6]; % mostly for electronic calib with varamp pulser
idx=find(T>TimeRange(1) & T<TimeRange(2));
T_=T(idx); A_=A(idx);


RawSigmaT=std(T_)
% all T_s are now inside time range
% clear T A

%------ cut lower amplitudes ---------------------------------------

As_=sort(A_);
[A_,T_]=cut(find(A_>As_(length(A_)*.02)),A_,T_);

% slice the amplitude spectrum and correct ---------------------------------------
NrBounds=13;	% 23 for many events
%NrBounds=23;	% 23 for many events

NrBins=NrBounds-1;
Sigmas_=zeros(1,NrBins);
Means=zeros(1,NrBins);
Centers=zeros(1,NrBins);

A__=[];T__=[];
%ACorr=[];TCorr=[];

[As_,I]=sort(A_);
Ts_=T_(I);
clear I
lAs_=length(As_);
coefs=[];

%Eff_TOF_Software=length(Ts_)/Events*100

Abnd=[linspace(1,lAs_/10*3/4,3) linspace(lAs_/10,lAs_,NrBounds-3)]; % subdivide first bin
%Abnd=[linspace(1,lAs_,NrBounds)]; 
for i=2:NrBounds
	Tt_=Ts_(Abnd(i-1)+1:Abnd(i));
	Aa_=As_(Abnd(i-1)+1:Abnd(i));
	MeansT_(i)=mean(Tt_);
	Centers(i)=mean(Aa_);  % centers of the A-bins
	if Centers(i) < 1000
% 		p=polyfit(Aa_,Tt_,2);
		p=polyfit(Aa_,Tt_,1);
	else
%		p=[0  0 polyfit(Aa_,Tt_,0)];
		p=[0 polyfit(Aa_,Tt_,0)];
	end
	Tt_=(Tt_-polyval(p,Aa_));

	coefs=[coefs; p];
%	[Mean, Sigma, Max]=stdfit(Tt_,'test');
	[Mean, Sigma, Max]=stdfit(Tt_);
	Sigmas(i)=Sigma;
	MeansT__(i)=Mean;
	A__=[A__; Aa_];
	T__=[T__; Tt_];

%figure(1);	plot(A_(idx),T_(idx),'.');	Sigmas(i)	Means(i)
%figure(2);	histf(T_(idx),50); pause;
end

A__=A__(20:length(A__));  % for some stupid reason this 10 points spoil the fit
T__=T__(20:length(T__));
Sigmas=sqrt(Sigmas.^2-SigmaStart^2);

TApp=mkpp(As_(Abnd(:))',coefs(:));
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
	axis([0 log10(2400) 0 100]); grid on;
else
	axis([0 2400 0 100]); grid on;
end
ylabel('rms tof');
xlabel('Amplitude');ylabel('sigma TOF');
set(gca,'UserData','rmsTOF');


figure
%subplot(2,1,1)
C=10.^(Centers)+min(Arpc)-1;    %A=log10(Arpc-min(Arpc)+1);
S=Sigmas./sqrt(2);
%S=Sigmas;
C=C*fCperDiv;
C=C/2;
plot(C,S,marker); hold on % convert to single chamber
axis([0 2400 0 max(S)]); grid on;
ylabel('Sigma tof/channel');
xlabel('Fast charge/channel (fC)');ylabel('sigma TOF');
set(gca,'UserData','rmsTOF');
logx; xaxis(10,2048*fCperDiv)

%subplot(2,1,2);
%histf(Arpc,0:10:2400);
%logx; xaxis(100,1000)

drawnow

CorrRmsT=sqrt(std(T__).^2-SigmaStart^2);

RawMeanT=mean(MeansT_)

% risetime plot ------------------------------------------------------------------
X=(MeansT_(2:NrBounds)-min(MeansT_(2:NrBounds)))/1000;
Y=min(Centers(2:NrBounds))./Centers(2:NrBounds);
%figure, plot(X,Y,'o',X,Y,'-');

% T-A plots ------------------------------------------------------------------
figure 

NrScat=rand(1,1000)*(length(T__)-1)+1;
subplot(2,2,1);
plot(A_(NrScat),T_(NrScat),'.'); hold on; ylabel('Tempo medido (ps)'); xlabel('log(Amplitude)')

plot(Centers,MeansT_,'or');
plot((Centers(:)*[1 1])',(MeansT_(:)*[1 1]+Sigmas(:)*[-1 1])','-r');
if length(findstr(titulo,'RPC'))
	if all(A<10)
		axis([0 log10(2400) mean(TimeRange)-1000 mean(TimeRange)+1000]);  % for module
	else
		axis([0 2400 TimeRange]);
	end
else
	if all(A<10)
		axis([0 log10(2400) mean(TimeRange)-1000 mean(TimeRange)+1000]);
	else
		axis([0 2400 TimeRange]);
	end

end

if all(A<10)
	plot(linspace(1,log10(2400),100),myppval(TApp,linspace(1,log10(2400),100)),'b'); % plot the interpolant
else
	plot(1:10:2400,myppval(TApp,1:10:2400),'b'); % plot the interpolant
end

title(titulo);
drawnow

subplot(2,2,2); 
plot(A__(NrScat),T__(NrScat),'.'); hold on; ylabel('Tempo corrigido (ps)');xlabel('log(Amplitude)')
plot(Centers,MeansT__,'or');
plot((Centers(:)*[1 1])',(MeansT__(:)*[1 1]+Sigmas(:)*[-1 1])','-r');
if length(findstr(titulo,'RPC'))
	if all(A<10)
		axis([0 log10(2400) -1000 1000]);  % for module
	else
		axis([0 2400 -1000 1000]);
	end
else
	if all(A<10)
		axis([0 log10(2400) -2000 2000]);
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

% Residual T -------------------------------------------------------
L=length(T__)/NrBins/5;
T=conv(T__,ones(1,L))/L;
T=T(L:length(T__));
%figure
%plot(T,'.'); hold on;
%plot(Abnd,MeansT__,'or');
%xlabel('amplitude (arb units)');
%ylabel('T(ps)');
%title([titulo ' - residual T (from sliding mean)' ]);
%ax=axis;
SigmaResidualT=std(T(length(T)/10:length(T)))
%text(min(ax(1:2)),max(ax(3:4)*8/10),['rms = ' num2str(SigmaResidualT) ' ps (excluding the origin)'])

%drawnow

% Residual Tstart correlation  -------------------------------------------------------

if ~exist('Varampdata') 

I=find(isnan(Tcorr));Tcorr(I)=0*I;

%disp('---- Residual Tstart correlation  -----');
%Tcorr=uncorr(Tcorr,Tstart);

%disp('---- Residual Amplitude correlation  -----');
%Tcorr=uncorr(Tcorr,Arpc);

%disp('---- Residual Baseline correlation  -----');
%Tcorr=uncorr(Tcorr,Ch22);

%disp('---- Residual scintilator amplitude correlation  -----');
%Tcorr=uncorr(Tcorr,Aap1);

end

% TOF sigma -------------------------------------------------------

[Mean, Sigma, Max]=stdfit(Tcorr,titulo);Sigma=Sigma
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

