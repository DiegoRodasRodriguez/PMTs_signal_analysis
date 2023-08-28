% DAQ events/spill

clear all
 

file = fopen('c:\DAQEvents.txt','w');
Data=[];

   filename=['dat_0167'];
   fullname=['D:\DATA.m\' filename '.mat'];
   
  %   return
   
   if exist(fullname)
      disp('===================================')
      load(fullname)
      Clock=U(:,29)/1e6;
      [InewSpill, LIspillNr]=findspills2(Clock,1);
      Data=[Data; [1, mean(diff(InewSpill))]];
  end
    
Data
fprintf(file,'%5.3f\t%5.3f\n', Data ); 
status = fclose('all');
%csvwrite('C:\DAQEvents.csv',a \t  b);
return


% X5-scan
Data=[];save collect Data % reset collect
for i=341:358	% só para
   filename=['dat_' sprintf('%04d',i)];
   if exist(['\gsi2003\data\' filename '.mat'])
      ana, collect
      mosaic, pause(4)
   end
end
csvwrite('X5-scan.csv',Data);
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


