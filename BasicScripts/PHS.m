% Pulse height spectrum


figure;
[N,X]=histf(Arpc,0:2400);
subplot(2,2,1); stairs(X,N); xaxis(0,2400);
title (titulo);ylabel('Counts');

[N,X]=histf([Arpc(:)' (0:2400)],0:2400);
subplot(2,2,2); stairs(X,N); logy;xaxis(0,2400); % for module

[N,X]=histf([Ch42(:)' (0:2400)],0:2400);
subplot(2,2,4); stairs(X,N); logy;xaxis(0,2400); % for module

[N,X]=histf(Ch42,0:2400);
subplot(2,2,3); stairs(X,N); xaxis(0,2400);
title (titulo);ylabel('Counts');

%figure;stairs(

drawnow;
return

figure
histf([Arpc; (1:max(Arpc))'],1:max(Arpc));
logx,logy,
drawnow


[N,X]=histf(Appc,100);
subplot(2,2,3); stairs(X,N);
%xaxis(0,1200);
xaxis(0,2400);
title ('PPC');xlabel('channels');ylabel('Counts');
subplot(2,2,4); stairs(X,N); logy;
%xaxis(0,1200);
xaxis(0,2400);
drawnow

