function getheader(fid)

if fid == (-1)
	error(['getheader: Could not open file: ' filename]);
end

for i=1:8 fgets(fid); end

return