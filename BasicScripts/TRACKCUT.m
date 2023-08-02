%A=ChargeX3(:,1);x=0:20:1000;
A=HitsX1(:,1);x=-10:10;
%A=PosX3(:,1);x=0:0.6:30;
%A=Aap;x=200:16:1000;
%A=Arpc;x=0:10:500;
%A=TCorr1;x=-1000:10:1000;


B=cut(find(Tcorr>400),A);
[NA,X]=histf(A,x);
[NB,X]=histf(B,x);

figure
%subplot(2,2,1);stairs(X,NA);hold on;stairs(X,NB);grid on
%subplot(2,2,3);stairs(X,NB./(NA+1));grid on

C=cut(find(Tcorr>-50&Tcorr<50),A);
[NC,X]=histf(C,x);

subplot(3,1,1);stairs(X,NA);hold on;stairs(X,NB);stairs(X,NC);grid on;
a=get(gca,'Children');set(a(2),'Color','blue');set(a(3),'Color','green');
yaxis(1,1000);logy,

%subplot(3,1,2);stairs(X,(NC+0)./(NB+0));grid on;title('central/tail')
subplot(3,1,2);stairs(X,(NB+0)./(NA+0));grid on;title('tail/all')
subplot(3,1,3);stairs(X,(NA)./(NC));grid on;title('all/central')
