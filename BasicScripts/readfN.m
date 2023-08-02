function readf(varargin)
%  Read data from file, detecting the number of colummns
%
% varargin 1 = filename
% varargin 2 = pathIn
% varargin 3 (if exist) = pathOut
 
filename=varargin(1);
pathIn=varargin(2);

if(nargin==3)
    pathOut=varargin(3);
else
    pathOut=pathIn;
end



fid = fopen([pathIn{1} filename{1} '.dat'],'r');
if fid == (-1)
	error(['readf: Could not open file: ' filename{1}]);
end

D=diff(find(fscanf(fid, '%c',1000)==' '));
fclose(fid);

I=find(D==13);
N=length(find((D(1:(I(1)))~=1)));




fid = fopen([pathIn{1} filename{1} '.dat'],'r');
if fid == (-1)
	error(['readf: Could not open file: ' filename{1}]);
end

U=fscanf(fid, '%g ',[N,inf])';U=U(:,1:N-1);

fclose(fid);

disp(['=== file ' filename{1} ' with ' num2str(N-1) ' columns and ' num2str(length(U(:,1))) ' rows']);


eval(['save ' pathOut{1} filename{1}]);



return

