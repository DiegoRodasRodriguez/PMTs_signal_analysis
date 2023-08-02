function printw(filename,fig)
% function printw(filename [,fig])
% print 1 figure or all figures

disp('printw allways appends!');
disp('please confirm filename (touch any key or ctrl-C)!');
pause

if nargin==1
	fig=findobj('Resize','on');
end

fig=sort(fig);
for i=1:length(fig)
	u=sprintf('print -dps -append print\\%s -f%d',filename,fig(i))
	eval(u);
end
