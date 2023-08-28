% gain 

% must FINDSPILL before

Nbins=100;

%N=find(Ntriggers>25 & Ntriggers<35);
N=find(Ntriggers>1 & Ntriggers<4500);
N=N(1:(floor(length(N)/Nbins)*Nbins));
Appc=Appc(N);Arpc=Arpc(N);Tppc=Tppc(N);Trpc=Trpc(N); Tstart=Tstart(N);Tevent=Tevent(N);TppcT1=TppcT1(N);TrpcT1=TrpcT1(N); 
Ntriggers=Ntriggers(N);
if exist('AS1')
	AS1=AS1(N);
end

N=floor(length(N)/Nbins)*Nbins;
figure
subplot(2,3,1);
Ntr=mean(reshape(Ntriggers(1:N),N/Nbins,Nbins));
plot((1:Nbins)*N/Nbins,Ntr)
title('Number of triggers');
hold on
%plot(Ntriggers(),'.');
yaxis(0,max(yaxis));
drawnow

subplot(2,3,2);
hist(Ntriggers,50);

subplot(2,3,3);
Tt=mean(reshape(Tstart(1:N),N/Nbins,Nbins))*47;
plot((1:Nbins)*N/Nbins,Tt)
title('T0 - T1 time walk');
yaxis(1e4,max(yaxis));

subplot(2,3,4);
A=mean(reshape(Arpc(1:N),N/Nbins,Nbins));
plot((1:Nbins)*N/Nbins,A)
title('RPC gain walk');
hold on
%Arpc=Arpc(1:N)-interp1(0:99,A,(0:(N-1))/N*99,'linear')+mean(A);
%A=mean(reshape(Arpc(1:N),N/Nbins,Nbins));
%plot((1:Nbins)*N/Nbins,A)
yaxis(0,max(yaxis));
drawnow

% time

subplot(2,3,5);
T=mean(reshape(Trpc(1:N),N/Nbins,Nbins))*47;
plot((1:Nbins)*N/Nbins,T)
title('RPC time walk');
hold on
yaxis(0,max(yaxis));

subplot(2,3,6);
Ts=std(reshape(Trpc(1:N),N/Nbins,Nbins))*47;
plot((1:Nbins)*N/Nbins,Ts)
title('RPC rms jitter walk');
hold on
yaxis(0,max(yaxis));

%subplot(2,1,1);
%Trpc=Trpc(1:N)-interp1(0:99,T,(0:(N-1))/N*99,'linear')+mean(T);
%T=mean(reshape(Trpc(1:N),N/Nbins,Nbins));
%plot((1:Nbins)*N/Nbins,T)
%yaxis(0,max(yaxis));

%subplot(2,1,2);
%Ts=std(reshape(Trpc(1:N),N/Nbins,Nbins))*47;
%plot((1:Nbins)*N/Nbins,Ts)
%yaxis(0,max(yaxis));

drawnow


N=1:N;
Appc=Appc(N);Arpc=Arpc(N);Tppc=Tppc(N);Trpc=Trpc(N); Tstart=Tstart(N);Tevent=Tevent(N);TppcT1=TppcT1(N);TrpcT1=TrpcT1(N); 
Ntriggers=Ntriggers(N);
if exist('AS1')
	AS1=AS1(N);
end

