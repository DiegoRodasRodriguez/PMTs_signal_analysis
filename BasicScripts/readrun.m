function readrun(file_basename)
%read csv values from DAQ pc files
% auto checks number of fields

filename=['D:\camac\gsi2003\' file_basename '.dat']
%filename=['..\..\dat_file\' file_basename '.dat']

fid = fopen(filename,'r');
if fid == (-1)
   disp(['readrun: Could not open file: ' filename]);
   return
end

fgets(fid);fgets(fid); % first line is often blank

s=fgets(fid);
N=(length(s)-length(find(real(s)==32))-12)/4+1

if floor(N)~=N 
   error('Wrong N');
end

fclose(fid);

U=getcsv(filename,N); % for .dat BEWARE THAT NUMBER OF FIELDS MUST BE GIVEN

clear u N i;
I=find(file_basename == ' ');if(I>0);file_basename(I) = '+';end;

eval(['save D:\camac\gsi2003\' file_basename ]);

return
