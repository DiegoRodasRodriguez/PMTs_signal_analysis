function [varargout]=eraele(I,varargin);
%elimina los elementos de los vectores indexados por I
%              [variables de salida]=eraele(I,variables de entrada);
%
%
%
flag=ones(length(varargin{1}),1);
flag(I)=0;I=find(flag==1);
for j=1:(nargin)-1
varargout{j}=varargin{j}(I);
end
return

