function U=readrun_A(basedir, basename, ext, N)

filename = [basedir, '\', basename, '.', ext];

fid = fopen(filename,'r');
if fid == (-1)
   disp(['readrun: Could not open file: ' filename]);
   return
end

getheader(fid);
U=getcsv_A(fid,N)

clear N;
I=find(filename == ' ');if(I>0);filename(I) = '_';end;
%If blank space in the name then replace by '_'

eval(['save ' basedir, '\', basename, '.', ext '.mat']);

return
