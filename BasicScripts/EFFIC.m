% eficiency vs. cut


%idx=find(Tppc>max(Trpc)-3);TofHardCutPpc=max(Appc(idx));
%idx=find(Trpc>max(Trpc)-3);TofHardCutRpc=max(Arpc(idx));
%pedRange=10:max([TofHardCutRpc TofHardCutPpc])+20; % the channels range to display in fig 2

pedRange=10:80; % the channels range to display in fig 2

figure; 
[Nrpc,Xrpc]=histf(Arpc,1:max(Arpc));
EFFrpc=(1-cumsum(Nrpc)/Events)*100;
plot(Xrpc(pedRange),EFFrpc(pedRange)); grid on
title ('--PPC, __ RPC');xlabel('amplitude cut(channels)');ylabel('% Counts above cut');
drawnow

[Nppc,Xppc]=histf(Appc,1:max(Appc));
EFFppc=(1-cumsum(Nppc)/Events)*100;
hold on
plot(Xppc(pedRange),EFFppc(pedRange),'--'); grid on
drawnow

for i=0:3
	index=max(find(Xppc<OFFSETppc+(i*SIGMAppc)));
	plot(Xppc(index),EFFppc(index),'*');
	index=max(find(Xrpc<OFFSETrpc+(i*SIGMArpc)));
	plot(Xrpc(index),EFFrpc(index),'+');
end

%index=max(find(Xppc<TofHardCutPpc));plot(Xppc(index),EFFppc(index),'*');
%index=max(find(Xppc<TofHardCutRpc));plot(Xrpc(index),EFFrpc(index),'o');

yaxis(90,100)
drawnow
