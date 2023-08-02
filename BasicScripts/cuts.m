function [varargout]=cut(c,I,varargin);
%  Get the data indexed by I, if c='get'
%  Del the data indexed by I, if c='del'
%  if a variable is a matrix then rows are indexed.


%in the future only c='del'
if(c=='del' | c==1)
J=1:length(varargin{1});
J(I)=0;
I=find(J);
end

for j=1:(nargin)-2
s=size(varargin{j});
   if(s(1)>1 & s(2)>1)
      varargout{j}=varargin{j}(I,:);
      %disp(['=== Warnig the variable ' num2str(j) ' is a Matrix']);
      %disp(['=== check is you want to index rows.']);
   else
   varargout{j}=varargin{j}(I);
   end
end

return
