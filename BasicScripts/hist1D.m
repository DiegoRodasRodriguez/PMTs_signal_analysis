%HIST1D
%Fast algorithm for ploting 1D histograms of given range x with range giving 
%the bin edges.

function [N,x]=hist1D(y,x)
N = histc(y, x);

if(all(size(x)~=size(N))), N = N'; end
if nargout == 0,      stairs(x,N); end
end