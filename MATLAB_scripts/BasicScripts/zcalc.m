
d=0.05;     %glass width 0.5 mm
h=0.086;    %PCB height (Valencia)
N=8;        %8 gaps
g=0.025;    %0.25 mm
eps_r   = 7;  %epsilon for glass
eps_pcb = 5;
c=30;       %cm  

eps_o     = 8.85e-5; %nF/cm

Z = 50;     %ohms

width = 1/Z*sqrt( ((N/2+1)*d*eps_pcb + N/2*g*eps_r*eps_pcb +2*h*eps_r)*((N/2+1)*d + N/2*g + 2*h) )/ sqrt(eps_r*eps_pcb) /(eps_o*c)
