% brute force

BeginOfSpill=[1;find(diff(Tevent)<0)+1];
Ntriggers=Tevent*0;
for i=1:(length(BeginOfSpill)-1)
	r=BeginOfSpill(i):BeginOfSpill(i+1);
	Ntriggers(r)=r*0+length(r)-1;
end

return

% older runs

figure
%histf(SpillID,1:max(SpillID));ylabel('events/spill');xlabel('# Spill')

% brute force
N=length(SpillID);
ID=zeros(1,N);
k=1;
for i=2:N
	if SpillID(i-1)==SpillID(i)
		ID(i)=k;
		k=k+1;
	else
		k=1;
	end

end
% there is a peak at ID=1. Dont consider
M=min(max(ID),80);
%M=80; % cut some strange things
MeanArpc=zeros(1,M);
MeanAppc=zeros(1,M);

for i=2:M
	MeanArpc(i-1)=mean(Arpc(find(ID==i)));
	MeanAppc(i-1)=mean(Appc(find(ID==i)));
end
plot(MeanArpc); hold on;
plot(MeanAppc,'--');
ylabel('Mean signal (channels)'); xlabel('# event in spill (No spills in ALICE)');
text(3, 500, 'RPC ___');
text(3, 450, 'PPC ---');
drawnow;
% cut beginning of spill

%idx=find(ID>30 & ID<80);
idx=find(ID<15);


Arpc=Arpc(idx);
Trpc=Trpc(idx);
Appc=Appc(idx);
Tppc=Tppc(idx);

disp(sprintf('Cuts made - %d events after cut',length(Arpc)));

return

U=Arpc(1:2^min(log2(length(Arpc)),16));
f=fft(U-mean(U));
ff=abs(fftshift(f));
ff=conv(ff,ones(1,100));
N=length(ff);
subplot(2,1,1)
plot(ff(N/2+1:N/2+3000));
title('RPC');

U=Appc(1:2^min(log2(length(Appc)),16));
f=fft(U-mean(U));
ff=abs(fftshift(f));
ff=conv(ff,ones(1,100));
N=length(ff);
subplot(2,1,2)
plot(ff(N/2+1:N/2+3000));
title('PPC');
