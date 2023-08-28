function M=getcsv_D(fid,Nfields)
% The fields must be separated by comma.
% FIXME: More generic reader will be provided!

if fid == (-1)
	error('getcsv_D: Could not open file!');
end

M=fscanf(fid, '%g %*c',[Nfields,inf]);

fclose(fid);
return
