function dev = genTvsDx_withZZ(param)
% generate one piezo-electric zigzag switch device
% Component: lense fiber -> fixed coupler -> coupler held by zigzag with
% reflector at the end
% WTJ & FM, 181001
% TODO: generate the correct release window, necessary for electrodes


%%
w_end = param.w_end;
w_mid = param.w_mid;
zz_w = param.zz_w;
g_mid = param.g_mid;   % gap between middle couplers
L2 = param.L2;      % typically 10 um
L = param.L;        % typically 10 um
g = param.g;        % gap between zigzag arms
dx = param.offset; % offset of both waveguides wrt each other in x direction
g_metal = param.g_metal;
%%
dev = gpack.Group(0,0,{});
layer_LN = 'M1_Neg';
layer_holes = 'M1_Pos_holes';

P_tether.L = 3.5;
P_tether.w_tether = 1;
P_tether.w2 = 2.4;         % waveguide width in connection with tether

L_tetherEnd2Edge = 10;
P_coupler.w = 0.8;          % waveguide width


%% generate reclector
reflector=genReflector();
reflector.translate([-dx;0]);
dev.addelement(reflector);

d_tetherCenter2ZZEnd = 3.5;
y_tetherCenter =  - L_tetherEnd2Edge - P_tether.L;
r_tetherArm = Rect(0,y_tetherCenter, d_tetherCenter2ZZEnd*2, P_tether.w_tether,...
    'base','center');
r_tetherArm.translate([-dx;0]);
r_tetherArm.layer = layer_LN;
dev.addelement(r_tetherArm);

%% generate pairs of zigzags
zz = genSingleZZBender_181002(struct('w',zz_w,'g', g ,'L2', L2,'w_joint',0.1,'L',L));
zz.translate([-dx;y_tetherCenter]);
dev.addelement(zz); %pair near reflector


%%
w_end = P_coupler.w_end;
[cc, ys] = genCouplerCoupler_noTether();
cc.translate([0; g_mid]);
dev.addelement(cc);

%% generate small electrodes on anchor
%extract ports
port_RL=zz.elements{1}.elements{2}.port1+[zz.x+zz.elements{1}.x;zz.y+zz.elements{1}.y];
port_RR=zz.elements{1}.elements{3}.port1+[zz.x+zz.elements{1}.x;zz.y+zz.elements{1}.y];
port_LL=[zz.elements{2}.elements{2}.port2(1,:);zz.elements{2}.elements{2}.port1(2,:)]+[zz.elements{2}.x;zz.elements{2}.y];
port_LR=[zz.elements{2}.elements{3}.port2(1,:);zz.elements{2}.elements{3}.port1(2,:)]+[zz.elements{2}.x;zz.elements{2}.y];

%right zigzag electrodes
param_elec.dy=6;
param_elec.w_pad=4;
param_elec.zz_w=zz_w;
param_elec.port_L=port_RL;
param_elec.port_R=port_RR;
[small_R_elec,port_B_RL,port_B_RR]=genSmallAnchorElectrodes(param_elec);
%small_R_elec.translate([-2*dx;0]);
dev.addelement(small_R_elec);

%left zigzag electrodes
param_elec.port_L=port_RL; %need to mirror electrodes later
param_elec.port_R=port_RR;
[small_L_elec,port_B_LL,port_B_LR]=genSmallAnchorElectrodes(param_elec);
small_L_elec.mirror([0;0],[0;1]);
small_L_elec.translate([-2*dx;0]);
dev.addelement(small_L_elec);

%% generate second pair of zigzag for symmetry
zz2 = genSingleZZBender_181002(struct('w',zz_w,'g', g ,'L2', L2,'w_joint',0.1,'L',L));
%zz2.rotate(180);
zz2.translate([0;g_mid+L_tetherEnd2Edge+P_tether.L]);
dev.addelement(zz2); %pair near edge coupler
%tether arm
d_tetherCenter2ZZEnd = 3.5;
y_tetherCenter =  g_mid+L_tetherEnd2Edge + P_tether.L;
l_tetherArm = Rect(0,y_tetherCenter, d_tetherCenter2ZZEnd*2, P_tether.w_tether,...
    'base','center');
l_tetherArm.layer = layer_LN;
dev.addelement(l_tetherArm);


%% make end coupler at (0,0)
dev.translate(0,-ys(end)-g_mid);

%% Add bonding pads
w_pad = 200;
w_pad_big = 260;
y_pad1 = -200-w_pad/2;
d_pad = 525;
y_pad2 = y_pad1-d_pad;
y_pad3 = y_pad2-d_pad;
pad1 = genBondPad(0,y_pad1,w_pad);
pad2 = genBondPad(0,y_pad2,w_pad);
pad3 = genBondPad(0,y_pad3,w_pad_big);
dev.addelement(gpack.Group(0,0,{pad1,pad2,pad3}));



%% Add big wires
w_wire=5;
d_wire2pad=10;
thin_feed=Rect([0;(y_pad2+y_pad3)/2+w_wire/2],w_pad_big,w_wire);
thin_feed.layer='metal_BB';
dev.addelement(thin_feed);

%right wires
v_endEL=port_B_RL;
v_endER=port_B_RR;
v1 = v_endEL; v2 = [0;y_pad1]; x1 = v1(1);
wire_ZZ2EL = Wire({v1,v2},w_wire);
% wire from ZZ2ER to line1
x2 = w_pad/2 + d_wire2pad;
v1 = v_endER; v2= [v_endER(1)+5;v_endER(2)];
v3=[v_endER(1)+5;y_pad1/2];
v4 = [x2; y_pad1/2];
v5 = [x2; (y_pad2+y_pad3)/2]; %v3 = [x2; (y_pad1+y_pad2)/2];
%v4 = [-w_pad/2 - d_wire2pad, (y_pad1+y_pad2)/2];
%v5 = [-w_pad/2 - d_wire2pad, y_line1];
wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);%wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);
g_wires = gpack.Group(0,0,{wire_ZZ2EL, wire_ZZ2ER});%wire_ZZ1EL, wire_ZZ1ER, rect_line1, rect_line2
g_wires.layer = 'metal_BB';
g_BB = gpack.Group(0,0,{g_wires});
dev.addelement(g_BB);

%left wires
v_endEL=port_B_LL;
v_endER=port_B_LR;
%wire from ZZ2EL to pad1
v1 = v_endEL+[2*dx;0]; v2 = [v_endEL(1)+2*dx;y_pad1/2]; x1 = v1(1);
x2 = w_pad/2 + d_wire2pad;
v3 = [x2;y_pad1/2];
v4 = [x2;(y_pad1+y_pad2)/2];
v5 = [-x2-w_wire/2;(y_pad1+y_pad2)/2];
wire_ZZ2EL = Wire({v1,v2,v3,v4,v5},w_wire);
% wire from ZZ2ER to line1
x2 = max(x1+ d_wire2pad,w_pad/2 + d_wire2pad);
x3=x2+10;
v1 = v_endER+[2*dx;0]; v2 = [v_endER(1)+2*dx+4; v1(2)];
v3 = [v_endER(1)+2*dx+4;y_pad1/4];
v4 = [x3;y_pad1/4];
v5 = [x3; y_pad2]; %v3 = [x2; (y_pad1+y_pad2)/2];
v6 = [0; y_pad2];
%v5 = [-w_pad/2 - d_wire2pad, y_line1];
wire_ZZ2ER = Wire({v1,v2,v3,v4,v5,v6},w_wire);%wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);

g_wires = gpack.Group(0,0,{wire_ZZ2EL, wire_ZZ2ER});%wire_ZZ1EL, wire_ZZ1ER
g_wires.layer = 'metal_BB';
g_wires.mirror([0;0],[0;1]);
g_BBL = gpack.Group(0,0,{g_wires});
dev.addelement(g_BBL);

%% add undercut mask
edge_mask = Rect([0;10],w_pad_big,20);
edge_mask.layer = 'M2_undercutMsk';
edge_mask.translate(0,ys(end)+g_mid);
dev.addelement(edge_mask);
%% add label
[ind_pos, ind_neg] = genStringPolygon(sprintf('A%d',param.ind),...
        0.4);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
ind.rotate(-90);
ind.translate(0,-50);

dev.addelement(ind);

%% functions for components
    function dev = genReflector()
    % generate reflector with end located at (0,0)
    % 

    dev = gpack.Group(0,0,{});

    % mirror cell from FBH GA 180707
    P_mirror = struct('a', 557.6e-3, 'w', 1083e-3, 'amp', 81.9e-3, 'isFBH',true,...
        'hx', 330e-3, 'hy', 859e-3);

    L_mirror2tether = 3.5;   % end of mirror cell region (roughly) to starting point of the tether


    P_coupler.L_taper = 3.5;     % taper length from NB width to WG width
    P_coupler.w = 0.8;          % waveguide width
    P_coupler.N_remove = 0;
    P_coupler.N_mirror = 18;
    P_coupler.N_coupler = 8;
    P_coupler.isfullcav = 0;
    P_coupler.l_extra = 0;     % extra coupler length into the buffer region
    P_coupler.maker = @FBHMaker_F;

    P_coupler.P_defect = P_mirror;
    P_coupler.P_mirror = P_mirror;
    P_coupler.N_trans = 0;
    P_coupler.holeLayer = layer_holes;
    P_coupler.slabLayer = layer_LN;
    l_tot = P_coupler.N_mirror * P_coupler.P_mirror.a;
    % apply correction
    P_coupler.w_end = w_end;
    P_coupler = correct_LNNB16_old(P_coupler);
    % waveguide with width P_coupler.w startes after the
    % P_coupler.L_taper transition region and ends at the edge between
    % exposure box 1 and box 2 (box 2 = buffer box for dicing)

    % generate coupler Polygon        
    % assuming no displacement or rotation of nanobeam
    %         dy = (l_NB_unscaled/2 - P_mirror.a/2);
    dy = 0;
    % right-bottom corner of NB before 90 deg rotation:
    x0 = (l_tot - P_mirror.a/2) - dy;         
    y0 = (P_mirror.w - 2 * P_mirror.amp)/2;

    % the following Ls and Ws include the tether for the coupler
    Ls = [P_coupler.L_taper, L_mirror2tether ,...
        P_tether.L, P_tether.L, L_tetherEnd2Edge];
    Ws = [P_coupler.w, P_coupler.w,P_tether.w2,P_coupler.w, w_mid];
    P_coupler.Ls = [];
    P_coupler.Ws = Ws;
    N_totalcoords = 2*(length(Ws) + 1);
    coordCell = {};
    coordCell{1} = [x0;y0];
    coordCell{N_totalcoords} = [x0;-y0];
    for ii = 1:length(Ws)
        x_tmp = x0 + sum([Ls(1:ii)]);
        coordCell{ii+1} = [x_tmp;Ws(ii)/2];
        coordCell{N_totalcoords-ii} = [x_tmp;-Ws(ii)/2];
    end
    % generate coupler slab to the right end 
    % this is done here to prevent zooming the NBs changing the coupler
    % parameters
    couplerSlab = Polygon(coordCell,P_coupler.slabLayer); 
    couplerSlab.rotate(90); 
    dev.addelement(couplerSlab); 

    coupler = nanobeamcoupler(rmfield(P_coupler,{'L_taper','w','l_extra','w_end'}));
    coupler.rotate(90);         
    coupler.translate(0, - dy );

    dev.addelement(coupler);
    
    % generate writefield box
    w_wf = 62;
    rect_wf = Rect(0,w_wf/2 - P_mirror.a/2,w_wf,w_wf,'base','center');
    rect_wf.layer = 'M1_Field';
    dev.addelement(rect_wf);
    
    % generate curing region
    w_cure = 6;
    h_cure = 66;
    rect_cure = Rect(0, h_cure/2 - P_mirror.a/2, w_cure, h_cure,...
        'base','center');
    rect_cure.layer = 'cure';
    dev.addelement(rect_cure);
    
    w_mask=10;
    optics_mask = Rect([0;h_cure/2 - P_mirror.a/2-4],w_mask,h_cure);
    optics_mask.layer = 'M2_undercutMsk';
    dev.addelement(optics_mask);
    
    % now the whole structure is vertical and the lowest points are at y = 0
    % going to translate the structure so that the coupler end is at (0,0)
    l_total = x0 + sum(Ls);
    dev.translate([0;-l_total]);
    end

    function [dev, ys] = genCouplerCoupler_noTether()
        % generate coupler-to-coupler structure
        
        Ws = [w_mid,P_coupler.w,P_tether.w2,P_coupler.w, w_end];
        Ls = [L_tetherEnd2Edge,P_tether.L, P_tether.L, L_tetherEnd2Edge];
        ys = cumsum(Ls);
        v1 = [Ws(1)/2;0]; v2 = [Ws(2)/2; ys(1)]; v3 = [Ws(3)/2;ys(2)]; v4 = [Ws(4)/2;ys(3)]; v5 = [Ws(5)/2;ys(4)];
        nx = @(v)([-v(1);v(2)]);
        p_couplerCoupler = Polygon({v1,v2,v3,v4,v5,nx(v5),nx(v4),nx(v3),nx(v2),nx(v1)});
        
        L_tether = 25;  % separation of anchors
        %rect_tether = Rect(0,ys(2),L_tether,P_tether.w_tether,'base','center');
        r = 15/2;
        %c_anchor_r = Circ(L_tether/2, ys(2), r);
        %c_anchor_l = Circ(-L_tether/2, ys(2), r);
        devHD = gpack.Group(0,0,{p_couplerCoupler});
        devHD.layer = layer_LN;
        %devLD = gpack.Group(0,0,{c_anchor_l, c_anchor_r});
        devLD.layer = 'M1_LD';
        dev = gpack.Group(0,0,{devHD});
    end


    function [dev, ys] = genCouplerCoupler()
        % generate coupler-to-coupler structure
        
        Ws = [w_mid,P_coupler.w,P_tether.w2,P_coupler.w, w_end];
        Ls = [L_tetherEnd2Edge,P_tether.L, P_tether.L, L_tetherEnd2Edge];
        ys = cumsum(Ls);
        v1 = [Ws(1)/2;0]; v2 = [Ws(2)/2; ys(1)]; v3 = [Ws(3)/2;ys(2)]; v4 = [Ws(4)/2;ys(3)]; v5 = [Ws(5)/2;ys(4)];
        nx = @(v)([-v(1);v(2)]);
        p_couplerCoupler = Polygon({v1,v2,v3,v4,v5,nx(v5),nx(v4),nx(v3),nx(v2),nx(v1)});
        
        L_tether = 25;  % separation of anchors
        rect_tether = Rect(0,ys(2),L_tether,P_tether.w_tether,'base','center');
        r = 15/2;
        c_anchor_r = Circ(L_tether/2, ys(2), r);
        c_anchor_l = Circ(-L_tether/2, ys(2), r);
        devHD = gpack.Group(0,0,{p_couplerCoupler, rect_tether});
        devHD.layer = layer_LN;
        devLD = gpack.Group(0,0,{c_anchor_l, c_anchor_r});
        devLD.layer = 'M1_LD';
        dev = gpack.Group(0,0,{devHD, devLD});
    end

end
