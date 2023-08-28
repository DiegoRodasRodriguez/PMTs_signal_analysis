% DAQ events/spill

clear all
Data=[];

eff=1;
pos_res=1*0;
if eff

file = fopen('c:\Eff.txt','w'); 

%167:211 == RPC I
%for ii=[223:247 285:359]	% para cint. a 14 cm
%for ii=[167:219]	  % RPC I   [167:218]

%for ii=[220:255]      % RPC III
    
%for ii=[339:359]       % RPC V    

%for ii=[219:3]

%for ii=[224:327]

for ii=[228:246]


%for ii=[219:330]


%for ii=[223:224]
% for ii=[339:358]   
    if (ii>=253 & ii<=275)
    continue;
    end    
    
    if ( ii==168 )
    continue;
    end    
    
    filename=['dat_' sprintf('%04d',ii)];
    fullname=['D:\DATA.m\' filename '.mat'];

   if exist(fullname) 
       disp(num2str(ii))
       load(fullname)
       Clock=U(:,29)/1e6;
       [InewSpill, LIspillNr]=findspills2(Clock,1); %pause(2);
       ana;
%      return
    
       hold off; 
       spills; pause(2);devrate
  %     Data=[Data; [ii, round(mean(diff(InewSpill))), Nout, Nlost, total(1), Neff, Ndouble, Ncross, devrate, meansigma, stdviation, Eff_after_external_cuts]];
       Data=[Data; [ii, round(mean(diff(InewSpill))), Nout, Nlost, total(1), Neff, Ndouble, Ncross, devrate, Doubles, Eff_after_external_cuts]];
       close all;     
  end
   
end

Data
fprintf(file,'%6d \t %6d \t %6d \t %6d \t %6d \t %6d \t %6d \t %6d \t %4.3f \t %4.3f \t %4.3f \n', Data' ); 

status = fclose('all');
%csvwrite('DAQEvents.csv',Data);

return

end

if pos_res

    file = fopen('c:\res.txt','w'); 
      
    %for ii=[167:211]   
    %for ii=[219:330] 
    %for ii=[228:243] 
    for ii=[339:358] 
        
    if (ii>=253 & ii<=275)
    continue;
    end    
    
    if ( ii==168 )
    continue;
    end    
    
    filename=['dat_' sprintf('%04d',ii)];
    fullname=['D:\DATA.m\' filename '.mat'];

   if exist(fullname) 
       disp(num2str(ii))
       load(fullname)
%       Clock=U(:,29)/1e6;
%       [InewSpill, LIspillNr]=findspills2(Clock,1); %pause(2);
       ana;
       mosaic; pause(2);
%      return
       hold off; 
%       spills; pause(2);devrate
       Data=[Data; [ii, res_x, mean_x]];
       close all;     
  end
   
end

Data
fprintf(file,'%6d \t %4.3f \t %4.3f \n', Data' ); 

status = fclose('all');
%csvwrite('DAQEvents.csv',Data);
    
return

end
    
    
% X5-scan
Data=[];%save collect Data % reset collect


filename=['dat_' sprintf('%04d',i)];
fullname=['D:\DATA.m\' filename '.mat'];

for ii=167:219	% só para
    
filename=['dat_' sprintf('%04d',ii)];
fullname=['D:\DATA.m\' filename '.mat'];

    if exist(fullname)
      ana, %collect
      mosaic, pause(4)
      Data=[Data; [ii, Sigmafinal]];   
  end
end


fprintf(file_s1,'%6d\t%6.3f\n', Data' ); 
status = fclose('all');

%csvwrite('X5-scan.csv',Data);

return


% X3-scan
Data=[];save collect Data % reset collect
for i=228:243	% só para
   filename=['dat_' sprintf('%04d',i)];
   if exist(['\gsi2003\data\' filename '.mat'])
      ana, collect
      mosaic, pause(4)
   end
end
csvwrite('X3-scan.csv',Data);
return

% batch dat->mat (precisa mudar directorio tambem em readrun)
for i=1:359	% só para
   filename=['dat_' sprintf('%04d',i)];
   if exist(['C:\Paulo\Hades\BEAM2003\' filename '.dat'])
      readrun(filename);    
   end
end
return


