clear all
close all

% runs 103 and 104
%x=650:50:1000;
%zx=[627 1300 2295 3000 2844 1982 996 330];
%y=900:25:1150;
%zy=[226 404 2727 3208 3304 3300 3267 3053 2600 565 201];
%peak_rate=3000/4; % Hz/cm2

% runs 105 and 106 107 % no peak!
x=600:100:1200;
zx=[1349 1440 1459 1501 1707 1700 1651];
y=500:100:1200;
zy=[1885 1852 1746 1658 1628 1482 1379 900];
peak_rate=800/4; % Hz/cm2


x=(x-mean(x))/10;
y=(y-mean(y))/10;
zx=zx/max(zx);
zy=zy/max(zy);


xi=linspace(min(x),max(x),100);
zix=interp1(x,zx*peak_rate,xi,'spline');

yi=linspace(min(y),max(y),100);
ziy=interp1(y,zy*peak_rate,yi,'spline');


plot(xi,zix,'-',yi,ziy,'-');
hold on;
plot(x,zx*peak_rate,'o',y,zy*peak_rate,'o');

figure
theta=linspace(0,2*pi,100);
%polar(theta,r);

Ax=[-20 -20 20 20 -20];Ay=[-20 20 20 -20 -20];
fill(Ax*1.05,Ay*1.05,'r'); hold on
fill(Ax,Ay,'w');
xaxis(-30,30);yaxis(-30,30);
axis('image')

area=zeros(1,10);

for i=9:-1:3
xx=[min(find(zix>i/10)) max(find(zix>i/10))];
yy=[min(find(ziy>i/10)) max(find(ziy>i/10))];
X0=mean(xx); Y0=mean(yy);
if i==9
	X00=X0;Y00=Y0;
end
a=(xx(2)-X0)/2; b=(yy(2)-Y0)/2;
r=a*b./sqrt(a^2+(b^2-a^2)*cos(theta).^2);
X=cos(theta).*r+X0-X00; 
Y=sin(theta).*r+Y0-Y00;
plot(X,Y); hold on

area(i)=a*b*pi/100; % cm2
end

set(2,'Position',get(2,'Position')*1.5)

figure
i=9:-1:3;
dArea=[area(9) diff(area(i))];
bar(i/10*peak_rate,dArea);
xlabel('Rate (Hz/cm2)');ylabel('Iluminated area (cm2)');
