close all
clear all
filename='dat_0417';
load 'data\dat_0418';
Q=U(:,21);
%Q=U(:,1);


load 'data\dat_0358';  % pedestal run
P=U(:,21);

figure;%subplot(1,2,1);
[N,X]=histf(Q,0:10:2100);stairs(X,N);
hold on;
[Np,Xp]=histf(P,0:10:2100);stairs(Xp,Np);
set(max(get(gca,'children')),'color',[0 0 1])
logy
xlabel('Q')
%title(filename);
%subplot(1,2,2);
figure
plot(X,(1-cumsum(N)/sum(N))*100);
hold on
plot(Xp,(1-cumsum(Np)/sum(Np))*100,'b');
grid
xlabel('Q')
