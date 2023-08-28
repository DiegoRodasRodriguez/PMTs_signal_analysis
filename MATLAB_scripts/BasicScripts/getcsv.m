function M=getcsv(filename,Nfields)
% The fields must be all SPACE delimited. CR is ignored.

%
% open the file 
%
fid = fopen(filename,'r');
if fid == (-1)
	error(['getcsv: Could not open file: ' filename]);
end

M=fscanf(fid, '%g ',[Nfields,inf])';

fclose(fid);
return

% Read delimited format 
%
eol = 10;         % End Of Line char
loc = [1 1];      % starting location of return matrix
line = fgets(fid); % get the 1st line, if any...
%
% read till eof
%
while(line ~= [ -1 ])
	disp (line)
	% get next line of file
	line = fgets(fid); 
end
% close file
fclose(fid);