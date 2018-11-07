function dev = genTvsDx_noZZ_three_devs(param)
% generate T vs dx device with no zigzags
% Component: lense fiber -> fixed coupler -> coupler held by zigzag with
% reflector at the end
% WTJ & FM, 181105
% TODO: generate the correct release window, necessary for electrodes


%%
w_end = param.w_end;
w_mid = param.w_mid;
g_mid = param.g_mid;   % gap between middle couplers
dx = param.offset; % offset of both waveguides wrt each other in x direction
d = 10; %distance between devices
%%
dev = gpack.Group(0,0,{});
layer_LN = 'M1_Neg';
layer_holes = 'M1_Pos_holes';

P_tether.L = 3.5;
P_tether.w_tether = 1;
P_tether.w2 = 2.4;         % waveguide width in connection with tether

L_tetherEnd2Edge = 10;
P_coupler.w = 0.8;          % waveguide width


%% generate reflector
reflector=genReflector();
reflector.translate([-dx;0]);
dev.addelement(reflector);

%%
w_end = P_coupler.w_end;
[cc, ys] = genCouplerCoupler_noTether();
cc.translate([0; g_mid]);
dev.addelement(cc);

%% copy and translate to left and right
dev_left = dev.copy();
dev_right = dev_left.copy();

% tmp_dev_left=dev_left.elements{1}.elements;
% dev_left.elements{1}.elements=[];
% dev_left.elements{1}.elements{1}=tmp_dev_left{1};
% dev_left.elements{1}.elements{2}=tmp_dev_left{2};
% dev_left.elements{1}.elements{3}=tmp_dev_left{4};


dev_left.translate([-d;0]);
dev_right.translate([d;0]);
dev.addelement(dev_left);
dev.addelement(dev_right);

%% generate lower long tether
d_tetherCenter2ZZEnd1 = 2*d;
y_tetherCenter1 =  - L_tetherEnd2Edge - P_tether.L;
r_tetherArm = Rect(0,y_tetherCenter1, d_tetherCenter2ZZEnd1*2, P_tether.w_tether,...
    'base','center');
r_tetherArm.translate([-dx;0]);
r_tetherArm.layer = layer_LN;
dev.addelement(r_tetherArm);
%% add upper long tether
d_tetherCenter2ZZEnd = 2*d;
y_tetherCenter =  g_mid+L_tetherEnd2Edge + P_tether.L;
l_tetherArm = Rect(0,y_tetherCenter, d_tetherCenter2ZZEnd*2, P_tether.w_tether,...
    'base','center');
l_tetherArm.layer = layer_LN;
dev.addelement(l_tetherArm);

%% generate writefield box
w_wf = 62.5;
rect_wf = Rect(0,-3,w_wf,w_wf,'base','center');
rect_wf.layer = 'M1_Field';
dev.addelement(rect_wf);
%% make end coupler at (0,0)
dev.translate(0,-ys(end)-g_mid);

%% generate upper anchors
%right
y0 = y_tetherCenter;%zz_LN.port1(2,1);
x1 = d_tetherCenter2ZZEnd-0.5*P_tether.w_tether-0.25;%zz_LN.port1(1,1);
x2 = d_tetherCenter2ZZEnd+0.5*P_tether.w_tether-0.25;%zz_LN.port1(1,2);
dy = 6;
x3 = x2 + dy;
x4 = x1 - dy;
y1 = y0 - dy;
y2 = y1 - dy;
p_tether = Polygon({[x2;y0],[x3;y1],[x3;y2],[x4;y2],[x4;y1],[x1;y0]});
p_tether.layer = 'M1_Neg';
r_rotCenter = [(x1+x2)/2;y0];
%p_tether.translate(0, zzp.w/2/sqrt(2));
p_tether.rotate(r_rotCenter,90);
dev.addelement(p_tether);
% add extra disk for tether/anchor
r = 15/2;
c_anchor = Circ((x1+x2)/2,y0 - 6-3, r,'M1_LD');
c_anchor.rotate(r_rotCenter,90);
dev.addelement(c_anchor);

%left
y0 = y_tetherCenter;%zz_LN.port1(2,1);
x1 = -d_tetherCenter2ZZEnd-0.5*P_tether.w_tether+0.25;%zz_LN.port1(1,1);
x2 = -d_tetherCenter2ZZEnd+0.5*P_tether.w_tether+0.25;%zz_LN.port1(1,2);
dy = 6;
x3 = x2 + dy;
x4 = x1 - dy;
y1 = y0 - dy;
y2 = y1 - dy;
p_tether = Polygon({[x2;y0],[x3;y1],[x3;y2],[x4;y2],[x4;y1],[x1;y0]});
p_tether.layer = 'M1_Neg';
r_rotCenter = [(x1+x2)/2;y0];
%p_tether.translate(0, zzp.w/2/sqrt(2));
p_tether.rotate(r_rotCenter,-90);
dev.addelement(p_tether);
% add extra disk for tether/anchor
r = 15/2;
c_anchor = Circ((x1+x2)/2,y0 - 6-3, r,'M1_LD');
c_anchor.rotate(r_rotCenter,-90);
dev.addelement(c_anchor);

%% generate lower anchors
%right
y0 = y_tetherCenter1;%zz_LN.port1(2,1);
x1 = d_tetherCenter2ZZEnd1-0.5*P_tether.w_tether-0.25;%zz_LN.port1(1,1);
x2 = d_tetherCenter2ZZEnd1+0.5*P_tether.w_tether-0.25;%zz_LN.port1(1,2);
dy = 6;
x3 = x2 + dy;
x4 = x1 - dy;
y1 = y0 - dy;
y2 = y1 - dy;
p_tether = Polygon({[x2;y0],[x3;y1],[x3;y2],[x4;y2],[x4;y1],[x1;y0]});
p_tether.layer = 'M1_Neg';
r_rotCenter = [(x1+x2)/2;y0];
%p_tether.translate(0, zzp.w/2/sqrt(2));
p_tether.rotate(r_rotCenter,90);
p_tether.translate([-dx;0]);
dev.addelement(p_tether);
% add extra disk for tether/anchor
r = 15/2;
c_anchor = Circ((x1+x2)/2,y0 - 6-3, r,'M1_LD');
c_anchor.rotate(r_rotCenter,90);
c_anchor.translate([-dx;0]);
dev.addelement(c_anchor);

%left
y0 = y_tetherCenter1;%zz_LN.port1(2,1);
x1 = -d_tetherCenter2ZZEnd1-0.5*P_tether.w_tether+0.25;%zz_LN.port1(1,1);
x2 = -d_tetherCenter2ZZEnd1+0.5*P_tether.w_tether+0.25;%zz_LN.port1(1,2);
dy = 6;
x3 = x2 + dy;
x4 = x1 - dy;
y1 = y0 - dy;
y2 = y1 - dy;
p_tether = Polygon({[x2;y0],[x3;y1],[x3;y2],[x4;y2],[x4;y1],[x1;y0]});
p_tether.layer = 'M1_Neg';
r_rotCenter = [(x1+x2)/2;y0];
%p_tether.translate(0, zzp.w/2/sqrt(2));
p_tether.rotate(r_rotCenter,-90);
p_tether.translate([-dx;0]);
dev.addelement(p_tether);
% add extra disk for tether/anchor
r = 15/2;
c_anchor = Circ((x1+x2)/2,y0 - 6-3, r,'M1_LD');
c_anchor.rotate(r_rotCenter,-90);
c_anchor.translate([-dx;0]);
dev.addelement(c_anchor);
%% add label
[ind_pos, ind_neg] = genStringPolygon(sprintf('P%d',param.ind),...
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
    
%     % generate writefield box
%     w_wf = 62;
%     rect_wf = Rect(0,w_wf/2 - P_mirror.a/2,w_wf,w_wf,'base','center');
%     rect_wf.layer = 'M1_Field';
%     dev.addelement(rect_wf);
    
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

    function [dev, ys] = genCouplerCoupler_noTether()
        % generate coupler-to-coupler structure
        
        Ws = [w_mid,P_coupler.w,P_tether.w2,P_coupler.w, w_end];
        Ls = [L_tetherEnd2Edge,P_tether.L, P_tether.L, L_tetherEnd2Edge];
        ys = cumsum(Ls);
        v1 = [Ws(1)/2;0]; v2 = [Ws(2)/2; ys(1)]; v3 = [Ws(3)/2;ys(2)]; v4 = [Ws(4)/2;ys(3)]; v5 = [Ws(5)/2;ys(4)];
        nx = @(v)([-v(1);v(2)]);
        p_couplerCoupler = Polygon({v1,v2,v3,v4,v5,nx(v5),nx(v4),nx(v3),nx(v2),nx(v1)});
        
        devHD = gpack.Group(0,0,{p_couplerCoupler});
        devHD.layer = layer_LN;
        dev = gpack.Group(0,0,{devHD});
    end



end
