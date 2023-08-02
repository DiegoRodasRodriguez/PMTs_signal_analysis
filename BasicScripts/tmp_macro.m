close all;
clear all;

segmentations=[1 2 3 4]

delta_t=0.5;

D=30*2/3*delta_t/40.0*segmentations;

Pinterf=(1-D).^2 

plot(segmentations,Pinterf,'s'); hold on;

delta_t=0.25;

D=30*2/3*delta_t/40.0*segmentations;

Pinterf=(1-D).^2 

plot(segmentations,Pinterf,'o');

delta_t=0.1;

D=30*2/3*delta_t/40.0*segmentations;

Pinterf=(1-D).^2 

plot(segmentations,Pinterf,'^');

figure;

voleff=0.72./(segmentations*2/3)*tan(pi/180*20);

plot(segmentations, voleff,'^');