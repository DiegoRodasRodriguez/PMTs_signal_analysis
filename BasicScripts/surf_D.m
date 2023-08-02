%Not-yet perfect approximation to a 2D histogram filled by using
%bin values and contents
%DGD (19/01/14)
%
%1) Color code in 2D is ok, but values are stored by MATLAB at vertexes :(.
%Values can be grafically retrieved by clicking on the low-left corner of
%the corresponding bin.
%
%2) White lines appear in the 2D representation when maximum and minimum
%values differ significantly. Can be solved by adjusting minimum and maximum 
%values and generally avoiding zeros. Possibly re-setting the color code
%here will also work.
%
%The input are the values at the bin centers and the corresponding mesh

function surf_D(Xmesh, Ymesh, N_xy)

% surf(Xmesh, Ymesh, N_xy, 'FaceColor', 'flat');
% view(2); colorbar;
% figure;
s_xy = size(N_xy);
l_y  = s_xy(1);
l_x  = s_xy(2);

N_xy_  = [zeros(1,l_x); N_xy;  zeros(1,l_x)];
N_xy__ = [zeros(l_y+2,1), N_xy_, zeros(l_y+2,1)];
X      = Xmesh(1,:);
pitchX = X(length(X))-X(length(X)-1);
Y      = Ymesh(:,1);
pitchY = Y(length(Y))-Y(length(Y)-1);

%Smart way, cause the sign of pitchY accounts for cases where the order of X, Y is descendent
[Xmesh_p, Ymesh_p] = meshgrid([X(1)-pitchX, X, X(length(X))+pitchX], [Y(1)-pitchY; Y; Y(length(Y))+pitchY]);

surf(Xmesh_p-pitchX/2, Ymesh_p-pitchY/2, N_xy__, 'FaceColor', 'flat');

minX = min(X) - abs(pitchX/2);
maxX = max(X) + abs(pitchX/2);
minY = min(Y) - abs(pitchY/2);
maxY = max(Y) + abs(pitchY/2);
view(2); colorbar;
xaxis(minX,maxX); yaxis(minY,maxY);