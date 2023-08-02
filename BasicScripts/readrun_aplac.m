function U1=readrun_aplac(basedir, basename, ext, N)
%read csv values from DAQ pc files
% auto checks number of fields

filename = [basedir, '\', basename, '.', ext];

fid = fopen(filename,'r');
if fid == (-1)
   disp(['readrun: Could not open file: ' filename]);
   return
end

%Eleven lines is the APLAC header
for i=1:11
fgets(fid); % first line is often blank
end

if(N==8) 
Time = fscanf(fid, '%*s %f %*s %*s %*s %*s %*s %*s %*s %*s %*s',[1, inf])'

%start again
fclose(fid);
fid  = fopen(filename,'r');
for i=1:11
fgets(fid); % first line is often blank
end
Volt = fscanf(fid, '%*s %*s %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s %f %*s',[N,inf])';

end


if(N==4) 
Time = fscanf(fid, '%*s %f %*s %*s %*s %*s %*s',[1, inf])'

%start again
fclose(fid);
fid  = fopen(filename,'r');
for i=1:11
fgets(fid); % first line is often blank
end
Volt = fscanf(fid, '%*s %*s %*s %f %*s %f %*s %f %*s %f %*s',[N,inf])';

end


if(N==3) 
Time = fscanf(fid, '%*s %f %*s %*s %*s %*s',[1, inf])'

%start again
fclose(fid);
fid  = fopen(filename,'r');
for i=1:11
fgets(fid); % first line is often blank
end
Volt = fscanf(fid, '%*s %*s %*s %f %*s %f %*s %f %*s',[N,inf])';

end

if(N==2) 
Time = fscanf(fid, '%*s %f %*s %*s %*s',[1, inf])'

%start again
fclose(fid);
fid  = fopen(filename,'r');
for i=1:11
fgets(fid); % first line is often blank
end
Volt = fscanf(fid, '%*s %*s %*s %f %*s %f %*s',[N,inf])';

end

if(N==1) 
Time = fscanf(fid, '%*s %f %*s %*s',[1, inf])'

%start again
fclose(fid);
fid  = fopen(filename,'r');
for i=1:11
fgets(fid); % first line is often blank
end
Volt = fscanf(fid, '%*s %*s %*s %f %*s',[N,inf])';

end



U1=[Time Volt];

fclose(fid);

eval(['save ' basedir, '\', basename ' U1']);

return
