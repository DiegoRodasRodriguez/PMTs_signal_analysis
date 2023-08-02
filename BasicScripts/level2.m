clear all
close all

load tmp2

tr=tr-mean(tr);
tp=tp-mean(tp);

%plot(tr,tp,'.');

t=150;

tp=tp+t;
t02=[];t04=[];t06=[];

figure
A=[-10;3];
for k=1:length(tr);
	p=polyfit([1 2],[tr(k) tp(k)],1);
	t=polyval(p,A);
	t04(k)=polyval(p,-4);
	t02(k)=polyval(p,-2);
	t06(k)=polyval(p,-6);

if k<100
	plot(t,A);hold on
	plot([tr(k) tp(k)],[1 2],'o');
end
end


hfitg(t02,40,1);
hfitg(t04,40,1);
hfitg(t06,40,1);
