% monte carlo to simulate the combination of independent detectors not fully efficient

clear all
close all

disp('-----------------------------')

Eff=0.75
N=2
Events=10000
data=zeros(1,Events);

Ineff=1-Eff;

for i=1:Events
   fired=0;t=0;
   for j=1:N
      if rand<Eff, fired=fired+1;t=t+randn(1,1);end
   end
   if fired, data(i)=t/fired; else data(i)=-20; end
end

data=data*100;
GlobalEff=length(find(data>-2000))/Events*100

histf(data,-1000:10:1000); 
hfitg(data,100,1,-1000,1000);yaxis(1,5000); logy
stdfit(data,' ');
