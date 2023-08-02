%Under development

function stairs2D(XY,N)

%        Nhits_mesh(length(Xmm_mesh)-i+1,j) = Nhits;
% UFF: does not look like the right way
%
%

X=XY(1,:);  %FIXME: this or the other dimension?. Better check for
%a rectangle.
Y=XY(:,1);
%Note: selection of a given vertex gives bin-content up-right.

%Sad necessity. Do faster?
N_ = zeros(size(XY));
for i=1:length(X), N_(i,:) = N(length(X)-i+1,:); end

h=bar3(X, N_, 1);

for i = 1:length(h)
    zdata = get(h(i),'ZData');
    set(h(i),'CData',zdata);
    set(h,'EdgeColor','k'); 
end

view(2);
colorbar;
set(gca, 'XTickLabel',  X); %How to do this?
%Easiest possibly to define a figure with the right axes,
%Then get it and set it. Nasty workaround once again