function U=readrun_D(basedir, basename, ext, N)
%read csv values from DAQ pc files
% auto checks number of fields

filename = [basedir, '\', basename, '.', ext];

fid = fopen(filename,'r');
if fid == (-1)
   disp(['readrun: Could not open file: ' filename]);
   return
end

%fgets(fid);fgets(fid); % first line is often blank
fclose(fid);

%s=fgets(fid);
%N=ceil((length(s)-length(find(real(s)==32))-12)/4+1);
%N=6
%disp('number of columns must be fixed!')

%Fixme!. This part I do not understand.

%if floor(N)~=N 
%   error('Wrong N');
%end

U=getcsv(filename,N); % for .dat BEWARE THAT NUMBER OF FIELDS MUST BE GIVEN
%csvread

%clear u N i; 
%I=find(filename == ' ');if(I>0);filename(I) = '+';end;

eval(['save ', basedir, '\', basename, ' U']);

return
