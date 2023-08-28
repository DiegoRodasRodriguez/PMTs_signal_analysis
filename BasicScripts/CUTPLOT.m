% disp(' ------- gen corrs ---------------------');
Aap=(Ch13+Ch14)/2;
Ar1p=(Ch15+Ch16)/2;
Ar2p=(Ch17+Ch18)/2;

Aam=(Ch13-Ch14)./(Aap+1);
Ar1m=(Ch15-Ch16)./(Ar1p+1);
Ar2m=(Ch17-Ch18)./(Ar2p+1);

Tap=(Ch25+Ch26)/2;
Tr1p=(Ch27+Ch28)/2-Tap;
Tr2p=(Ch29+Ch30)/2-Tap;

Tam=(Ch25-Ch26);
Tr1m=(Ch27-Ch28);
Tr2m=(Ch29-Ch30);

if ~exist('Tcorr'), return,end

figure
subplot(4,4,1);plot(Tcorr(1:1000),Aam(1:1000),'.');ylabel('Aam'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,2);plot(Tcorr(1:1000),Aap(1:1000),'.');ylabel('Aap'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,3);plot(Tcorr(1:1000),Tam(1:1000),'.');ylabel('Tam'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,4);plot(Tcorr(1:1000),Tap(1:1000),'.');ylabel('Tap'); xaxis(-1000,1000); zoom on; grid on

subplot(4,4,5);plot(Tcorr(1:1000),Ar1m(1:1000),'.');ylabel('Ar1m'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,6);plot(Tcorr(1:1000),Ar1p(1:1000),'.');ylabel('Ar1p'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,7);plot(Tcorr(1:1000),Tr1m(1:1000),'.');ylabel('Tr1m'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,8);plot(Tcorr(1:1000),Tr1p(1:1000),'.');ylabel('Tr1p-Tap'); xaxis(-1000,1000); zoom on; grid on

subplot(4,4,9);plot(Tcorr(1:1000),Ar2m(1:1000),'.');ylabel('Ar2m'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,10);plot(Tcorr(1:1000),Ar2p(1:1000),'.');ylabel('Ar2p'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,11);plot(Tcorr(1:1000),Tr2m(1:1000),'.');ylabel('Tr2m'); xaxis(-1000,1000); zoom on; grid on
subplot(4,4,12);plot(Tcorr(1:1000),Tr2p(1:1000),'.');ylabel('Tr2p-Tap'); xaxis(-1000,1000); zoom on; grid on

subplot(4,4,13);plot(Tcorr(1:1000),Ch20(1:1000),'.');ylabel('Astart'); xaxis(-1000,1000); zoom on; grid on
