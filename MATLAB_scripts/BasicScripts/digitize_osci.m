clear all
close all

disp('                                                             ')
disp('           *------------------------------------------------*')
disp('           *    Macro for digitizing signals from scope     *')
disp('           *------------------------------------------------*')
disp('                                                             ')
pause(1)

%HERE GIVE FILE NAME

%-------------------------analysis with scintillators
%base_name    = 'F:\Miguel_Data\Analysis_RPC\Jun08\';
%file_name    = '2scin_5cm_parallel_Na';              % 20 ns, 500 points
%file_name    = '2scin_5cm_parallel_Na_2.4kV_10Kevt_5ns'; %5 ns, 2500 points
%file_name    = '2scin_5cm_parallel_Na_2.4kV_10Kevt'; % 2 ns, 2500 points
%file_name    = '2scin_5cm_parallel_Na_2.2kV_2Kevt'; % 2 ns, 2500 points
%file_name    = '2scin_5cm_parallel_Na_2.5kV_2Kevt'; % 2 ns, 2500 points
%file_name    = 'Scint1Scint2_Par_Na22_2.4kV_080623_dat'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), 2.4kV

%%%%LAST
%file_name     = 'Scint1Scint2_Par_Na22_2.4kV_080626'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), th=-100 mV, 2.4 kV
%file_name     = 'Scint1Scint2_Par_Na22_2.5kV_080626'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), th=-100 mV, 2.5 kV
%file_name     = 'Scint1Scint2_Par_Na22_2.3kV_080626'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), th=-100 mV, 2.3 kV
%file_name     = 'Scint1Scint2_Par_Na22_2.2kV_080626'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), th=-100 mV, 2.2 kV
%file_name     = 'Scint1Scint2_Par_Na22_2.1kV_080626'; % 20 ns, 151 points. Real time, 4 channels (400 ps sampling), th=-100 mV, 2.1 kV
%file_name     = 'Scint2RPC_Na22_2.4kV3.0kV_080625'; % 20 ns, 151 points. Real time, 3+1 channels (400 ps sampling), th=-20 mV, 2.4 kV /3.0 kV
%file_name      = 'Scint1RPC_Na22_2.4kV3.0kV_080624'; % 20 ns, 151 points. Real time, 3+1 channels (400 ps sampling),  th=-20 mV, 2.1 kV /3.0 kV

%-------------------------analysis for RPC FEE
base_name     = 'F:\DANI_ELE_ANALISIS\';
%file_name     = 'Na10mV_200mV_081105';
%file_name     = 'CosmicRays_LogicTrigg_081106';
file_name     = 'CosmicRays_Avalanch_LogicTrigg_081113';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GIVE THE DIFFERENT Y-SCALE FACTORS FROM CHANNEL TO CHANNEL
%scintillators
%Vfact1        = 1; Vfact2        = 1; Vfact3        = 1; Vfact4        = 1;

%'CosmicRays_Avalanch_LogicTrigg_081113'; %avalanches (before ampli)
%Vfact1        = 1; Vfact2        = 1; Vfact3        = 100/500; Vfact4        = 1;
%'CosmicRays_LogicTrigg_081106'; %streamers (after ampli)
Vfact1        = 1; Vfact2        = 1; Vfact3        = 10/500; Vfact4        = 1;

%'Na10mV_200mV_081105';
%Vfact1        = 1; Vfact2        = 1; Vfact3        = 10/500; Vfact4        = 200/500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GIVE THE DIFFERENT IMPEDANCE IN OHMS
Z1            = 50; Z2            = 50; Z3            = 50; Z4            = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%GIVE THE NUMBER OF POINTS IN THE X-AXIS
npoints     =  500;  %Number of saved points per event

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%DEFINE THE TRIGGER MODE AND THRESHOLD

%trigger_mode  = 1;     %Two scintillators seen Na22 source
%trigger is defined as a a condition on all the channels
%Typical case: 4 channels for each PM going to the scope

trigger_mode  = 2;      %Trigger in one signal (LVDS or reference plastic) 
%trigger is defined as a condition in the trigger channel
%Typical case: 3 channels (TOT, RPC signal, trigger signal) going to the
%scope

threshold_1   = -20;   %Threshold in mV for the first detector. Sign does not matter
threshold_2   = -20;   %Threshold in mV for the second detector. Sign does not matter

%IMPORTANT NOTE: If the external trigger is set to 100 mV then this point ALWAYS
%appear in the wave function!!. The more convenient is then to set a common
%threshold for all the channels equal to the one used here. This should be
%more precise. Otherwise, an interpolation is needed!!!

%END OF INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%

%STRUCTURE OF OUTPUT FILE
%
%   T_A1          Q_A1          V_A1         TOT_A1     channel 1 
%   T_A2          Q_A2          V_A2         TOT_A2     channel 2
%   T_D1          Q_D1          V_D1         TOT_D1     channel 3
%   T_D2          Q_D2          V_D2         TOT_D2     channel 4
%
%   time at      charge in   amplitude    time over threshold at 50 %
%   threshold    window                   maximum
%
%
%

%FROM NOW ON NO NEED TO TOUCH!

file_name_out = [file_name,'_digi'];
extension     = '.mat';

file          = [base_name file_name extension];

eval(['load ', file]);

n_events    = length(U)/npoints -1  %Subtract the last one that sometimes is crap
                                    %FIXME: what if not integer??
%n_events    = 1000;
                                    
n_ineff=0;

V1(1:npoints)=0;
V2(1:npoints)=0;
V3(1:npoints)=0;
V4(1:npoints)=0;

for i=1:n_events

    i=i

    Xlow = (i-1)*npoints+1;
    Xup  = i*npoints;

    Time(1:npoints) = U(Xlow:Xup,1);              %in ps!
    V1(1:npoints)   = U(Xlow:Xup,2)*1000*Vfact1;  %in mV!
    V2(1:npoints)   = U(Xlow:Xup,3)*1000*Vfact2;
    V3(1:npoints)   = U(Xlow:Xup,4)*1000*Vfact3;
    V4(1:npoints)   = U(Xlow:Xup,5)*1000*Vfact4;

    Time_range = (max(Time(1:npoints)) - min(Time(1:npoints)))/1000; %in ns

    %This corrects for the different scales, when in principle the scale
    %of the channel one is stored in the header file

    %Determine and subtract pedestals. Take the entries up to 2 points before the
    %trigger. FIXME: 10 points before the maximum is arbitrary. Check with
    %the analysis

    %     Pedestal determined by averaging values 10 points before the
    %     minimum
    %     Iped1 = min(find(V1(i,1:npoints)==min(V1(i,1:npoints))))-10;
    %     Iped2 = min(find(V2(i,1:npoints)==min(V2(i,1:npoints))))-10;
    %     Iped3 = min(find(V3(i,1:npoints)==min(V3(i,1:npoints))))-10;
    %     Iped4 = min(find(V4(i,1:npoints)==min(V4(i,1:npoints))))-10;

    %    if(Iped1<1) Iped1=3; disp('problem with pedestal! '); end
    %    if(Iped2<1) Iped2=3; disp('problem with pedestal! '); end
    %    if(Iped3<1) Iped3=3; disp('problem with pedestal! '); end
    %    if(Iped4<1) Iped4=3; disp('problem with pedestal! '); end

    %     V1(i,1:npoints)   = V1(i,1:npoints) - (sum(V1(i,1:Iped1-10)))/length(1:Iped1-10);
    %     V2(i,1:npoints)   = V2(i,1:npoints) - (sum(V2(i,1:Iped2-10)))/length(1:Iped2-10);
    %     V3(i,1:npoints)   = V3(i,1:npoints) - (sum(V3(i,1:Iped3-10)))/length(1:Iped3-10);
    %     V4(i,1:npoints)   = V4(i,1:npoints) - (sum(V4(i,1:Iped4-10)))/length(1:Iped4-10);

    V1(1:npoints)   = V1(1:npoints) - (sum(V1(1:20)))/length(1:20);
    V2(1:npoints)   = V2(1:npoints) - (sum(V2(1:20)))/length(1:20);
    V3(1:npoints)   = V3(1:npoints) - (sum(V3(1:20)))/length(1:20);
    V4(1:npoints)   = V4(1:npoints) - (sum(V4(1:20)))/length(1:20);

    if(trigger_mode==1 & (min(V1(1:npoints))>threshold_1 | min(V2(1:npoints))>threshold_1 | min(V3(1:npoints))>threshold_2 | min(V4(1:npoints))>threshold_2))
        n_ineff=n_ineff+1;
        Q_A1(i) = -1; Q_A2(i) = -1; Q_D1(i) = -1; Q_D2(i) = -1;
        A_A1(i) = -1; A_A2(i) = -1; A_D1(i) = -1; A_D2(i) = -1;
        T_A1(i) = -1; T_A2(i) = -1; T_D1(i) = -1; T_D2(i) = -1;
        continue;
    end

    if(trigger_mode==2 & (max(V1(1:npoints))<abs(threshold_1)))
        n_ineff = n_ineff+1;
        Q_A1(i) = -1; Q_A2(i) = -1; Q_D1(i) = -1; Q_D2(i) = -1;
        A_A1(i) = -1; A_A2(i) = -1; A_D1(i) = -1; A_D2(i) = -1;
        T_A1(i) = -1; T_A2(i) = -1; T_D1(i) = -1; T_D2(i) = -1;
        ToT_A1(i) = -1; ToT_A2(i) = -1; ToT_D1(i) = -1; ToT_D2(i) = -1;
        continue;
    end


    %FIXME: assumes that the minimum is corresponding to the
    %trigger. Probably true almost always.
    %Always in volts
    if(trigger_mode==1) A_A1(i)     = -(min(V1(1:npoints)));
    elseif(trigger_mode==2) A_A1(i) = (max(V1(1:npoints)));
    end
    A_A2(i)=-(min(V2(1:npoints)));
    A_D1(i)=-(min(V3(1:npoints)));
    A_D2(i)=-(min(V4(1:npoints)));

    %Fixme. This is equivalent to integrate by boxes in a.u. Get charge in
    %the right units and integrate with trapezoids at least

    % Integral by bars
    Q_A1(i)=-(sum(V1(1:npoints)))*Time_range/npoints/Z1; %In pC
    %     Q_A2(i)=-(sum(V2(i,1:npoints)))*Time_range/npoints/Z2;
    %     Q_D1(i)=-(sum(V3(i,1:npoints)))*Time_range/npoints/Z3;
    %     Q_D2(i)=-(sum(V4(i,1:npoints)))*Time_range/npoints/Z4;

    % Integral by trapezoids

    % V1_diff_tmp = diff(V1(i,1:npoints));
    % V2_diff_tmp = diff(V2(i,1:npoints));
    % V3_diff_tmp = diff(V3(i,1:npoints));
    % V4_diff_tmp = diff(V4(i,1:npoints));
    % V1_tmp      = V1(i,1:(npoints-1));
    % V2_tmp      = V2(i,1:(npoints-1));
    % V3_tmp      = V3(i,1:(npoints-1));
    % V4_tmp      = V4(i,1:(npoints-1));

    V_A1_int = @ (time_) interp1(Time(1:npoints)/1000, V1(1:npoints), time_);
    V_A2_int = @ (time_) interp1(Time(1:npoints)/1000, V2(1:npoints), time_);
    V_D1_int = @ (time_) interp1(Time(1:npoints)/1000, V3(1:npoints), time_);
    V_D2_int = @ (time_) interp1(Time(1:npoints)/1000, V4(1:npoints), time_);

    time_   = Time(1)/1000:Time_range/(npoints*10):Time_range; %Time in ns
    VA1_tmp  = V_A1_int(time_);
    VA2_tmp  = V_A2_int(time_);
    VD1_tmp  = V_D1_int(time_);
    VD2_tmp  = V_D2_int(time_);
    Q_A1(i)=-(sum(VA1_tmp))*Time_range/npoints/10/Z1; %In pC
    Q_A2(i)=-(sum(VA2_tmp))*Time_range/npoints/10/Z2;
    Q_D1(i)=-(sum(VD1_tmp))*Time_range/npoints/10/Z3;
    Q_D2(i)=-(sum(VD2_tmp))*Time_range/npoints/10/Z4;



    %   Q_A1(i)  = -quadl(V_A1_int, Time(i,2), Time(i,200))/Z1; %fC
    %   Q_A2(i)  = -quadl(V_A2_int, Time(i,2), Time(i,200))/Z2; %fC
    %   Q_D1(i)  = -quadl(V_D1_int, Time(i,2), Time(i,50))/Z3; %fC
    %   Q_D2(i)  = -quadl(V_D2_int, Time(i,2), Time(i,50))/Z4; %fC


    V1max  = max(max(VA1_tmp));

    %In ns
    ToT_A1(i) = time_(max(find(VA1_tmp>0.5*V1max))) -  time_(min(find(VA1_tmp>0.5*V1max)));
    ToT_A2(i) = -1; ToT_D1(i) = -1; ToT_D2(i) = -1;

    %Condition for inefficiency: at least one event below threshold. I
    %think this should not be done. It must be done in the analysis later
    %on. FIXME!

    %Note: This takes the first time the comparator fires and ignores the
    %rest!. No retriggers are considered

    Iup1 = min(find(V1(1:npoints)<threshold_1));    Ilow1  = Iup1-1;
    Iup2 = min(find(V2(1:npoints)<threshold_1));    Ilow2  = Iup2-1;
    Iup3 = min(find(V3(1:npoints)<threshold_2));    Ilow3  = Iup3-1;
    Iup4 = min(find(V4(1:npoints)<threshold_2));    Ilow4  = Iup4-1;

    %Fixme: condition

    if(Ilow1==0 | Ilow2==0 | Ilow3==0 | Ilow4==0)
        T_A1(i) =  0; T_A2(i) =  0; T_D1(i) =  0; T_D2(i) =  0;
        continue;
    end

    if(trigger_mode==1)
        T_A1(i) = interp1(V1(Ilow1:Iup1), Time(Ilow1:Iup1), threshold_1); T_A2(i) = interp1(V2(Ilow2:Iup2), Time(Ilow2:Iup2), threshold_1);
        T_D1(i) = interp1(V3(Ilow3:Iup3), Time(Ilow3:Iup3), threshold_2); T_D2(i) = interp1(V4(Ilow4:Iup4), Time(Ilow4:Iup4), threshold_2);
    elseif(trigger_mode==2)
        T_A1(i) = -1; T_A2(i) = -1;
        T_D1(i) = -1; T_D2(i) = -1;
    end
    %Fixme. This muxt be done better, with an spline probably

end

if(trigger_mode==1)

    Q2ADC_1 = npoints/3*abs(threshold_1)*10;    %This gives roughly the maximum charge
    Q2ADC_2 = npoints/3*abs(threshold_2)*10;    %This gives roughly the maximum charge

    Q_A1= 2048*Q_A1/Q2ADC_1; %FIXME: this goes in ADC units for simplifying next macro
    Q_A2= 2048*Q_A2/Q2ADC_1;
    Q_D1= 2048*Q_D1/Q2ADC_2;
    Q_D2= 2048*Q_D2/Q2ADC_2;

    Q_A1(find(Q_A1>=2048))=2048;
    Q_A2(find(Q_A2>=2048))=2048;
    Q_D1(find(Q_D1>=2048))=2048;
    Q_D2(find(Q_D2>=2048))=2048;

end
  
inefficiency = n_ineff/n_events

U=[];V1=[];V2=[];V3=[];V4=[];Time=[]; 

file=[base_name, file_name_out, extension];
eval(['save ', file]);

return;