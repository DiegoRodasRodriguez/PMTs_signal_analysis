%function to emulate the h/pl in paw
function [n1,x1]=hpl(data,npts,miny,maxy)
[n1,x1]=hist2(data,npts,miny,maxy);
stairs(x1,n1)
line([x1(1) x1(1)],[0 n1(1)])
