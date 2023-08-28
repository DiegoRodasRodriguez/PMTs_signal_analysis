% analise rate capability

clear all
close all

load 'data\rpc2139'; 
Arpc=Arpc1;

% cut some events
N=1:25000; rate=300 % Hz/cm2.
%N=25000:95000; 
%N=1:80000;
%N=100000:124000; rate=1200 % Hz/cm2.

% system parameters
rate=300 % Hz/cm2.
ADCbin=.005 % pC/bin

% chamber parameters (per gap)
d0=0.3	% mm
V0=5800/4	% V
plate=0.3/2	% cm consider only half of the glass thickness because this is per gap
Cplate=8.8e-14*4/plate 	% F/cm2 
Cgap=8.8e-14/plate*10	% F/cm2
Rplate=.1e12*plate 	% Ohm/cm2
TotalToFastRatio=10

% gas parameters Q=N0*G=exp(Ad+B*V)=C*exp(B*V)
C=7.31e-4 	 	% pC
B=4.35e-3  		% V-1



Arpc=Arpc(N); Tevent=Tevent(N);

clear Arpc1 Arpc2 Trpc1 Trpc2 Tstart

Events=length(Arpc)

Tevent_=(Tevent-min(Tevent))/10000; % convert to s

[BeamRate,Time]=histf(Tevent_,linspace(0,.450,50));
BeamRate=BeamRate/mean(BeamRate)*rate; % convert to rate

A=Time*0;
for i=1:(length(Time)-1)
i
	I=find(Tevent_>Time(i) & Tevent_<=Time(i+1));
	if length(I)
		A(i)=mean(Arpc(I));
	end
end
A=A*ADCbin*TotalToFastRatio/4; % convert to collected pC/gap

figure;
subplot(2,2,1);
stairs(Time,BeamRate);
xlabel('time in spill (s)'); ylabel('Counting rate (Hz/cm2)');

subplot(2,2,2);
stairs(Time,A);
xlabel('time in spill (s)');ylabel('charge/event/gap (pC)')
%yaxis(4,7);

subplot(2,2,3);
V=log(A/TotalToFastRatio/C)/B;
stairs(Time,V); 
yaxis(min(V),1600);
xlabel('time in spill (s)');ylabel('Voltage drop across plate (V)')

I=BeamRate*mean(A)*1e-12; % current in A
%I=I*0+1;
Impulse=[0 exp(-(Time)/Rplate/(Cplate+Cgap))];
Vv=conv(I,Impulse)/(Cplate+Cgap)*min(diff(Time));
hold on
plot(Time,V0-Vv(1:length(Time))+100,'b')
%plot(Time,V0-Vv(1:length(Time))+70,'b')

subplot(2,2,4);
plot(Time,I);
xlabel('time in spill (s)');ylabel('Current/cm2 (A)')
