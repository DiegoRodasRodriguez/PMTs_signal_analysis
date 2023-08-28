function out=colors(i)



c=['y' 'm' 'c' 'r' 'b' 'k'];
if(i>6)i=i-6*(floor((i-1)/6));end
out=c(i);
return

