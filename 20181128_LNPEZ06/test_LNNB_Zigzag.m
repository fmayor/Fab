function dev_all = test_LNNB_Zigzag(param)
% testing LNNB + zigzag for tunable kappa_e
% WTJ, 20181016

if nargin < 1
    param = struct('scale_NB', 1,...
        'd_nb2refl', 1.8,'d_nb2refl2', 2.0,'w_end', 0.26, 'tether_w2', [2.4],...
        'tether_L',[4],'tether_w_tether', [1],'isGenPS',true,...
        'zigzag_tether',0.05);
    param.ind = 0;
    param.isGenBB = true;
    param.g_metal = 0.15;
    param.w_zigzag = 0.45;
    param.g_zigzag = 0.3;
    param.L2 = 15;
end

w_end = param.w_end;
d_nb2refl = param.d_nb2refl;
d_nb2refl2 = param.d_nb2refl2;
scale_NB = param.scale_NB;
ind = param.ind;
tether_w2 = param.tether_w2;
tether_L = param.tether_L;
tether_w_tether = param.tether_w_tether;


% mirror cell from FBH GA 180707
P_mirror = struct('a', 557.6e-3, 'w', 1083e-3, 'amp', 81.9e-3, 'isFBH',true,...
    'hx', 330e-3, 'hy', 859e-3);
% best cost:
% P_defect = struct('a', 480.7e-3, 'w', 1344e-3, 'amp', 277.5e-3, 'isFBH',true,...
%     'hx', 352e-3, 'hy', 548.8e-3);
% ind = 207 from GA 180723_2
P_defect = struct('a', 450e-3, 'w', 1170e-3, 'amp', 39e-3, 'isFBH',true,...
    'hx', 334e-3, 'hy', 811e-3);

dev = nanobeam_DS_EC_ZZ2;
dev.P_defect = P_defect;
dev.P_mirror = P_mirror;
dev.N_trans = 12;
dev.N_mirror = 14; %17
dev.d_tether2edge = 10+7;
dev.P_coupler.l_extra = 0;    % coupler end at buffer box edge
dev.P_tether.L = 3.5;
dev.P_tether.w_tether = 1;
dev.d_nb2edge = 7;
dev.d_nb2refl = d_nb2refl;
dev.d_nb2refl2 = d_nb2refl2;
dev.w_buffer = 62.5;
dev.layer_pos = 'M1_Pos';
dev.layer_neg = 'M1_Neg';
dev.layer_wf = 'M1_Field';
dev.layer_holes = 'M1_Pos_holes';
dev.layer_anchor = 'M1_LD';
dev.w_end = w_end;
dev.P_undercutMsk.isGen = 1;
dev.P_undercutMsk.layer = 'M2_undercutMsk';
dev.P_undercutMsk.d = 4.5;
dev.P_tether.w2 = tether_w2;
dev.P_tether.L = tether_L;
dev.P_tether.w_tether = tether_w_tether;
dev.P_tether.w_tetherArm = 16.2;

dev.scale_NB = scale_NB;
dev.maker = @FBHMaker_F;
dev.P_coupler.maker = @FBHMaker_F;


dev.P_zigzag.layer = 'M1_MD';
if param.isGenBB
    dev.P_zigzag.layer_metalBB = 'metal_BB';
else
    dev.P_zigzag.layer_metalBB = 'trash';
end
dev.P_zigzag.w_tether = param.zigzag_tether;
dev.P_zigzag.g_metal = param.g_metal;
dev.P_zigzag.w = param.w_zigzag;
dev.P_zigzag.g = param.g_zigzag;
dev.P_zigzag.L2 = param.L2;

% apply correction
correct_LNNB16(dev);
% do not run, wait for adding bended edge coupler
dev.run(false);



%% add big beam circuits
if param.isGenBB
    layer_metalBB = 'metal_BB';
    w_pad = 200;
    h_pad = 300;
    % w_noLOR = 190;
    w_wire = 5;
    w_feeder = 8;
    w_total = 62.5*4;
    y_pad1 = -250-18.5-h_pad/2;
    d_pad = h_pad+130;
    d_wire2pad = (w_total-w_pad)/4;
    y_pad2 = y_pad1 - d_pad;
    
    w_pad3 = 350;
    
    y_pad3 = y_pad2 - d_pad;
    y_line1 = (y_pad2 + y_pad3)/2;
    y_line2 = y_pad3;

    p1 = genBondPad(0,y_pad1,w_pad,h_pad);
    p2 = genBondPad(0,y_pad2,w_pad,h_pad);
    p3 = genBondPad(0,y_pad3,w_pad3,h_pad);
    rect_line1 = Rect(0, y_line1, w_total*1.2, w_feeder,'base','center');
    rect_line2 = Rect(0,y_line2, w_total*1.2, w_wire,'base','center');

    % wires from device to pads
    % left electrodes connect to pad1 and pad2, right electrodes connect to line1
    % wire from ZZ2EL to pad1
    v1 = dev.vs_EL(:,2); v2 = [0;y_pad1]; x1 = v1(1);
    wire_ZZ2EL = Wire({v1,v2},w_wire);
    % wire from ZZ2ER to line1
    d12 = [70/sqrt(3);-70];
    x2 = max(x1+ d_wire2pad,w_pad/2 + d_wire2pad);
    v1 = dev.vs_ER(:,2); 
    v12 = v1 + d12;
    v2 = [x2; v12(2)];
    v3 = [x2; (y_pad1+y_pad2)/2]; v4 = [-w_pad/2 - d_wire2pad, (y_pad1+y_pad2)/2];
    v5 = [-w_pad/2 - d_wire2pad, y_line1];
    wire_ZZ2ER = Wire({v1,v12,v2,v3,v4,v5},w_wire);
    % wire from ZZ1EL to pad2
    x3 = x2 + d_wire2pad;
    v1 = dev.vs_EL(:,1); 
    v12 = v1 + d12;
    v2 = [x3; v12(2)]; v3 = [x3; y_pad2]; v4 = [0;y_pad2];
    wire_ZZ1EL = Wire({v1,v12,v2,v3,v4},w_wire);
    % wire from ZZ1ER to line1
    x4 = x3 + d_wire2pad;
    v1 = dev.vs_ER(:,1);
    v12 = v1 + d12;
    v2 = [x4; v12(2)]; v3 = [x4; y_line1];
    wire_ZZ1ER = Wire({v1,v12,v2,v3},w_wire);


    g_wires = gpack.Group(0,0,{rect_line1, rect_line2, wire_ZZ2EL, wire_ZZ2ER,...
        wire_ZZ1EL, wire_ZZ1ER});
    g_wires.layer = layer_metalBB;

    g_BB = gpack.Group(0,0,{p1,p2,p3, g_wires});

    %%
    dev_all = gpack.Group(0,0,{dev, g_BB});
else
    dev_all = gpack.Group(0,0,{dev});
end

g_ind = genLNNB_DS_EC_ind(ind);
dev_all.addelement(g_ind);

% dev_all.todxf('tmp');


%%

% function dev = genBondPad(x,y, w_bondPad)
% 
% if nargin < 3
%     w_bondPad = 190;
% end
%     w_bondPadWindow = w_bondPad - 10;
%     r_bondPad = Rect(0, 0, w_bondPad, w_bondPad, 'base','center');
%     r_bondPad.layer = layer_metalBB;
%     r_bondPadWin = Rect(0, 0, w_bondPadWindow, w_bondPadWindow, 'base','center');
%     r_bondPadWin.layer = 'M2_undercutMsk';
%     dev = gpack.Group(x,y,{r_bondPad,r_bondPadWin});
% end




function ind = genLNNB_DS_EC_ind(ii)

mkrPixelSize = 0.5;
[ind_pos, ind_neg] = genStringPolygon(sprintf('K%d',ii),...
        mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = -1;
dy = -75-1.5;
ind.rotate(-90);
ind.translate(dx,dy);

end




end
