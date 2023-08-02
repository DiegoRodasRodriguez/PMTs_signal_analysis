% tracking 
% uses the first cluster.

% positions in cm

z=200; % position of my chamber
X1Off=16.3; X2Off=14.3; X3Off=13.9;
Y1Off=12.8; Y2Off=12.3; Y3Off=12.4;
Z1=37.5; Z2=314; Z3=426;

PX1=PX1-X1Off;PX2=PX2-X2Off;PX3=PX3-X3Off;
PY1=PY1-Y1Off;PY2=PY2-Y2Off;PY3=PY3-Y3Off;

figure;subplot(2,3,1);hfitg(PX1,100,2),subplot(2,3,2);hfitg(PX2,100,2),subplot(2,3,3);hfitg(PX3,100,2),
subplot(2,3,4);hfitg(PY1,100,2),subplot(2,3,5);hfitg(PY2,100,2),subplot(2,3,6);hfitg(PY3,100,2),drawnow

k=(z-Z1)/(Z2-Z1);
x=PX1+(PX2-PX1)*k;
y=PY1+(PY2-PY1)*k;

x=x*4;y=y*4; % convert to mm

I=find(TrpcT1<2048);
figure;plot(x,y,'.r');hold on;plot(x(I),y(I),'.g');
axis([-50 50 -50 50]); 
drawnow

k=(Z2-Z1)/(Z3-Z1);
DistP2X= PX1+(PX3-PX1)*k-PX2;
DistP2Y= PY1+(PY3-PY1)*k-PY2;
DistP2=sqrt(DistP2Y.^2+DistP2Y.^2);
figure; subplot(2,2,1);histf(DistP2X,-2:.1:2);subplot(2,2,2);histf(DistP2Y,-2:.1:2);
subplot(2,2,3);histf(DistP2,0:.1:10);





