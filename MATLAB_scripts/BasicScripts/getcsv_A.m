function M=getcsv_A(fid,Nfields)
% The fields must be all SPACE delimited. CR is ignored.

if fid == (-1)
	error(['getcsv_A: Could not open file: ' filename]);
end

M=fscanf(fid, '%g ',[Nfields,inf])';

fclose(fid);
return
