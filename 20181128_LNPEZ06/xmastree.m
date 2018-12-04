%%
pk = imread('xmastree_181203.png');

BW_top = ~im2bw(pk,graythresh(pk)*0.5);
        
        
%%
charBoundaries = [];

[B, L, N, A] = bwboundaries(BW_top);

% figure;
% hold all;
% for ii = 1:length(B)
% % for ii = 1:N
%     tmp = B{ii};
%     plot(tmp(:,1),tmp(:,2));
% end

coord.pos = {B{1:N}};
if length(B) == N
    coord.neg = {};
else
    coord.neg = {B{N+1:end}};
end
charBoundaries = [charBoundaries coord];
save('coord_xmastree.mat', 'charBoundaries');
        
%%
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
g_pos.layer = 'POS';
g_neg.layer = 'NEG';
g_all = gpack.Group(0,0,{g_pos, g_neg});
% g_all.rotate(-90);
% g_all.zoom(4);
% g_all.todxf('pumpkin');


%% test bendpolygongroup
g_all.todxf('xmastree');
% bendpolygongroup(g_all, [-0;144], 160, [-1;0]);
% g_all.todxf('pumpkin_181102');
