%HIST3D
% 3D histogram for uniform binning
%
% Adapted from PF 29-03-13 by DGD (more streamlined)
% 
% Simple solution but not completelly elegant. Values are stored at vertexes
% (function surfc) and bin content is read-out at the up-right bin have.
% Also white lines created by usage of surf might be avoided in a nicer way
%

function N=hist3D(Mxy, xrange, yrange)
edges{1} = xrange;
edges{2} = yrange;

[Xmm_mesh_plot, Ymm_mesh_plot] = meshgrid(xrange, yrange);
  
N = hist3(Mxy, 'Edges', edges);
if nargout ==0
    surfc(Xmm_mesh_plot, Ymm_mesh_plot, N', 'FaceColor', 'flat');
    xaxis(min(xrange),max(xrange)); yaxis(min(yrange),max(yrange));
    %set(get(gca,'child'),'EdgeColor','none');
    %set(get(gca,'child'),'FaceColor','none');
    view(2);
    colorbar;
end