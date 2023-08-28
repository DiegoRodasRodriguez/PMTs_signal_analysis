function U=readrun_osci(basedir, basename, ext, data_type)
%read csv values from DAQ-osci files

if(data_type<4)
    filename         = [basedir, '\', basename, '.', ext];
    filename_header  = [basedir, '\', basename, '.h'];
elseif(data_type>=4)
    filename         = [basedir, '\', basename, '_dat', '.', ext];
    filename_header  = [basedir, '\', basename, '_hdr', '.', ext];
end

fid        = fopen(filename,'r');
fid_header = fopen(filename_header,'r')

%data_type=2;         %Little bug. Last entry of channel 4 overlaps with 1st entry of channel 1.
                      %This is corrected if data_type =2. Also the event
                      %header is of a different type
                     
%data_type=3;         %New header

%data_type=4;         %New header. Little bug belonging to data_type=2
                      %now corrected. Problem now is that first entry of
                      %channel 1 is very big!. Corrected

%THIS IS THE LATEST FORMAT                      
                      
%data_type=5;         %File with both scintillators                      
                      
%data_type=7;         %File with RPC and only 3 columns         


if fid == (-1)
   disp(['readrun: Could not open file: ' filename]);
   return
end

if fid_header == (-1)
   disp(['readrun: Could not open file header: ' filename_header]);
   return
end

%Find number of points
header = fscanf(fid_header, '%c', [1,inf]);

I=find(header==';');

if(data_type==2)
 %npoints is 7th field, separated by ;
 Inpoints_low = I(6)+1;
 Inpoints_up  = I(7)-1;
 % 
 npoints=str2num(header(Inpoints_low:Inpoints_up))
 
 %timebin is 10th field, separated by ;
 Itimebin_low = I(9)+1;
 Itimebin_up  = I(10)-1;
 %in ps 
 timebin=str2num(header(Itimebin_low:Itimebin_up))*1e+12
 
 %voltagebin is 14th field, separated by ;
 Ivoltbin_low = I(13)+1;
 Ivoltbin_up  = I(14)-1;
 %in mV. FIXME: is it
 voltbin=str2num(header(Ivoltbin_low:Ivoltbin_up))

elseif(data_type>2)
 %npoints is 7th field, separated by ;
 %Take the number before the comma.   
    
 Inpoints_low = I(6)+6;
 Inpoints_up  = I(7)-1;
 % 
 npoints=str2num(header(Inpoints_low:Inpoints_up))
 
 %timebin is 10th field, separated by ;
 Itimebin_low = I(9)+4;
 Itimebin_up  = I(10)-1;
 %in ps 
 timebin=str2num(header(Itimebin_low:Itimebin_up))*1e+12
 
 %voltagebin is 14th field, separated by ;
 Ivoltbin_low = I(13)+5;
 Ivoltbin_up  = I(14)-1;
 %in mV. FIXME: is it
 voltbin=str2num(header(Ivoltbin_low:Ivoltbin_up))  
end


 nchannel = 1;
 nevent   = 1;
 
 Ilow     = 1;
 Iup      = npoints;
 
 Ilow_M   = 1;
 Iup_M    = npoints;
 M=fscanf(fid, '%g %*c',[1, 8000*4*npoints]); %inf. Only 2000*4*npoints because of memory problems!
  
 while(Iup_M<length(M))
                    
     %Column for time!
     if (nchannel==1) U(Ilow:Iup,1) = (1:npoints)*timebin;   end

     %Column for voltages!
     U(Ilow:Iup,nchannel+1)    = M(Ilow_M:Iup_M)'*voltbin;
          
     if(data_type==7 & nchannel==1) 
         U(Ilow:Iup,3) = U(Ilow:Iup,2);  
         nchannel = nchannel +1;
     end    %Duplicate channel 1
     
     %Incremental step
     Ilow_M = Iup_M + 1;
     Iup_M  = Iup_M + npoints;

     if ( nchannel==4 )
         nevent   = nevent   + 1
         if(data_type==2 | data_type==3)
             U(Iup,nchannel+1) = U(Iup-1,nchannel+1);  %Patch the last entry of 4th channel
             Ilow_M = Ilow_M - 1;   
             Iup_M  = Iup_M  - 1;
             if(Iup_M<length(M)) M(Ilow_M) = M(Ilow_M+1); end %Patch M for channel 1 if next event is existing
         end
         
         %Bug, first entry in channel 1 is not ok
        if( (data_type==4 | data_type==5 | data_type==7) & Iup_M<length(M) ) M(Ilow_M) = M(Ilow_M+1); end
            
         Ilow     = Iup+1;
         Iup      = Iup + npoints;
         nchannel = 0;         
     end
          
     nchannel = nchannel + 1;
     
 end

nevent   = nevent  
 
fclose(fid);
fclose(fid_header);

clear u N i;
I=find(filename == ' ');if(I>0);filename(I) = '+';end;

eval(['save ' basedir, '\', basename,'.mat' ]);

return
