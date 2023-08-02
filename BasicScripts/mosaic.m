%Function for re-distributing figures over screen
%Diego Gonzalez-Diaz 07-01-2013
%Up and down offsets might slightly depend on the screen resolution. Here
%they have been optimized for a resolution of 1440x900.

function mosaic(Yoff_down,Yoff_up)

%Before Matlab 14
%windows = sort(get(0,'children'));

%Starting from Matlab 14
FigArray = get(0,'children');
for i=1:length(FigArray), NumArray(i) = sort(FigArray(i).Number); end  %#ok<AGROW>
windows = sort(NumArray);

screen  = get(0,'screenSize');

%Offset to avoid the windows bar
if nargin<1     
    Yoff_down = 20; 
end

Ncols   = ceil(sqrt(length(windows)));
Nrows   = ceil(length(windows)/Ncols);

%Running offset so that header of upper figures is always inside display
if nargin<2 
    Yoff_up = ceil(100/Nrows + 33*(Nrows-1)); 
end

Xsize   = floor((screen(3) - screen(1))/Ncols);
Ysize   = floor((screen(4) - screen(2) - (Yoff_down + Yoff_up))/Nrows);    %Define 40

for i=Nrows:-1:1
    for j=Ncols:-1:1
        Iwindow = (i-1)*Ncols + j;
        Xo = screen(1) + (j-1)*Xsize;
        Yo = screen(4) - i*Ysize + Yoff_down - Yoff_up;
        WindowPos = [Xo, Yo, Xsize, Ysize];
        if(Iwindow<=length(windows))
            set(windows(Iwindow),'visible','off');
            set(windows(Iwindow),'position', WindowPos);
            set(windows(Iwindow),'visible','on');
        end
    end
end

