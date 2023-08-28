function readf(filename,N)

%fid = fopen(['a:\' filename '.dat'],'r');
fid = fopen(['d:\' filename '.dat'],'r');
%%%fid = fopen(['c:\gsi2003\data\' filename '.dat'],'r');
%fid = fopen(['c:\dat_file\tofpet\' filename '.dat'],'r');
%fid = fopen(['c:\temp\' filename '.dat'],'r');
%fid = fopen(['c:\pastas~1\alberto\result~1\camara\' filename '.txt'],'r');
%fid = fopen(['c:\pastas~1\alberto\dat_file\' filename '.dat'],'r');
%fid = fopen(['c:\pastas~1\alberto\dat_file\camara\' filename '.dat'],'r');
%fid = fopen(['c:\pastas~1\alberto\dat_file\tofpet\vsenergy\' filename '.dat'],'r');
%fid = fopen(['c:\pastas~1\alberto\dat_file\gepid\' filename '.dat'],'r');
%fid = fopen(['c:\tmp\' filename],'r');



if fid == (-1)
	error(['readf: Could not open file: ' filename]);
end

U=fscanf(fid, '%g ',[N,inf])';

fclose(fid);


I=find(filename == ' ');filename(I)='+';
%eval(['save c:\data\' filename]);

%eval(['save h:\' filename]);
%eval(['save c:\pastas~1\alberto\dat_file\gepid\' filename]);
%eval(['save c:\dat_file\tofpet\' filename]);
eval(['save c:\gsi2003\data\' filename]);
%eval(['save c:\pastas~1\alberto\dat_file\tofpet\vsenergy\' filename]);
%eval(['save c:\pastas~1\alberto\dat_file\camara\' filename]);
%eval(['save c:\pastas~1\alberto\tesina\datos\matlab\' filename]);
%eval(['save c:\pastas~1\alberto\result~1\camara\' filename]);
%eval(['save e:\dat\' filename]);
%eval(['save a:\' filename]);

return

