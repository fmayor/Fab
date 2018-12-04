function [g_pos_out, g_neg_out] = genXmastree_1801203(w_tot)

if nargin < 1
    w_tot = 450;
end

load('coord_xmastree.mat');

coords = charBoundaries;
coord_pos = coords.pos;
coord_neg = coords.neg;
    d = 0;
    g_neg = gpack.Group(0,0,{});
        g_pos = gpack.Group(0,0,{});
for ii = 1:length(coord_pos)
    p = Polygon(coord_pos{ii});
    p.translate(0,d);
    g_pos.addelement(p);
end
for jj = 1:length(coord_neg)
    p = Polygon(coord_neg{jj});
    p.translate(0,d);
    g_neg.addelement(p);
end
% g_pos.layer = 'POS';
% g_neg.layer = 'NEG';
% g_all = gpack.Group(0,0,{g_pos, g_neg});
g_pos.rotate(-90);
g_pos.translate(-212,204);
g_pos.zoom(w_tot/450);
g_pos_out = gpack.Group(0,0,{g_pos});
g_neg.rotate(-90);
g_neg.translate(-212,204);
g_neg.zoom(w_tot/450);
g_neg_out = gpack.Group(0,0,{g_neg});


end