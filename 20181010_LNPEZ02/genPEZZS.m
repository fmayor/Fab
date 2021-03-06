function dev = genPEZZS(param)
% generate one piezo-electric zigzag switch device
% Component: lense fiber -> fixed coupler -> coupler held by zigzag with
% reflector at the end
% WTJ & FM, 181001
% TODO: generate the correct release window, necessary for electrodes


%%
w_end = param.w_end;
w_mid = param.w_mid;
g_mid = param.g_mid;   % gap between middle couplers
L2 = param.L2;      % typically 10 um
L = param.L;        % typically 10 um
g = param.g;        % gap between zigzag arms

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
dev.addelement(genReflector());

d_tetherCenter2ZZEnd = 3.5;
y_tetherCenter =  - L_tetherEnd2Edge - P_tether.L;
r_tetherArm = Rect(0,y_tetherCenter, d_tetherCenter2ZZEnd*2, P_tether.w_tether,...
    'base','center');
r_tetherArm.layer = layer_LN;
dev.addelement(r_tetherArm);

%% generate pairs of zigzags
zz = genSingleZZBender_181002(struct('w',0.3,'g', g ,'L2', L2,'w_joint',0.1,'L',L));
zz.translate([0;y_tetherCenter]);
dev.addelement(zz);

%%
w_end = P_coupler.w_end;
[cc, ys] = genCouplerCoupler();
cc.translate([0; g_mid]);
dev.addelement(cc);

%% make end coupler at (0,0)
dev.translate(0,-ys(end)-g_mid);

%% add label
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',param.ind),...
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
    P_coupler = correct_LNNB16(P_coupler);
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
    h_cure = 62;
    rect_cure = Rect(0, h_cure/2 - P_mirror.a/2, w_cure, h_cure,...
        'base','center');
    rect_cure.layer = 'cure';
    dev.addelement(rect_cure);
    
    % now the whole structure is vertical and the lowest points are at y = 0
    % going to translate the structure so that the coupler end is at (0,0)
    l_total = x0 + sum(Ls);
    dev.translate([0;-l_total]);
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
