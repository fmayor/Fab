classdef nanobeam_DS_EC_ZZ2 < gpack.Group
% nanobeam with edge coupler for lense fiber
% lense fiber -> coupler -> reflector, two nanobeams side-coupled to the
% reflector
% Vertically oriented
% WTJ, 20180523
% Add secondary mask for XeF2, change first mask to adapt to HSQ
% WTJ, 20180615
% 
% Modified from nanobeam_DS_EC for zigzag tunable-kappa LNNB
% Left NB not changed, right NB connected to two zigzags
% WTJ, 20181016
properties
    w_end       % coupler end width
    w_buffer    % buffer region width
    h_buffer    % buffer region depth
    P_tether    % tether parameters
    P_coupler   % coupler parameters
    P_undercutMsk   % undercut mask parameters
    P_zigzag    % zigzag parameters
    P_mirror
    P_defect
    N_trans
    N_mirror
    maker       % nanobeam unitcell maker
    d_nb2refl    % center distance between right nanobeam and reflector
    d_nb2refl2   % center distance between left nanobeam and reflector
    d_nb2edge   % distance from nanobeam to undercut box edge
    d_nb2tether % vertical distance from NB end to tether on coupler WG
    d_tether2edge   % vertical distance from tether to buffer region edge
    scale_NB    % scaling factor for nanobeams, also scales the mirrors of the reflector
    ind_NB1     % index of the nanobeam object
    ind_NB2
    ind_coupler
    l_NB_unscaled
    layer_pos   % layer name of area to expose (no resist, to etch)
    layer_holes
    layer_neg   % layer name of area not to expose
    layer_wf    % write field layer
    layer_cure  % layer for HSQ curing
    layer_anchor  % layer for anchoring, lower dose
    vs_EL       % output ports of left/positive electrode, ONLY CORRECRT WITHOUT FURTHER TRANSLATION/ROTATION
    vs_ER       % output ports of right/negative electrode
end
methods
    function obj = nanobeam_DS_EC_ZZ2(varargin)
        obj@gpack.Group(0,0,{});
        P_tether = obj.getDefaultTetherParams();
        P_coupler = obj.getDefaultCouplerParams();
        P_undercutMsk = obj.getDefaultUndercutMskParams();
        P_zz = obj.getDefaultZigzagParams();
        % from GA180413
        P_defect = struct('a',0.467,'w',0.926,...
            'hx',0.205,'hy',0.616);
        % mirror param from GA20180412
        P_mirror = struct('a',0.576,'w',0.926,...
        'hx',0.332,'hy',0.558);
        p = inputParser;
        padp = @(varargin)(p.addParameter(varargin{:}));
        padp('w_end',0.25, @isnumeric);
        padp('w_buffer',62.5, @isnumeric);      % one write field
        padp('h_buffer',15,@isnumeric);
        padp('d_nb2refl',2,@isnumeric);
        padp('d_nb2refl2',2,@isnumeric);
        padp('d_nb2edge',10,@isnumeric);
        padp('d_nb2tether',5,@isnumeric);
        padp('d_tether2edge',7,@isnumeric);
        padp('scale_NB',1,@isnumeric);
        padp('P_tether',P_tether,@isstruct);
        padp('P_coupler',P_coupler,@isstruct);
        padp('P_undercutMsk',P_undercutMsk,@isstruct);
        padp('P_zigzag',P_zz,@isstruct);
        padp('P_mirror',P_mirror,@isstruct);
        padp('P_defect',P_defect,@isstruct);
        padp('N_trans',16,@isnumeric);
        padp('N_mirror',18,@isnumeric);
        padp('maker',[]);
        padp('layer_pos','pos');
        padp('layer_holes','pos_holes');
        padp('layer_neg','neg');
        padp('layer_wf','field');
        padp('layer_cure','cure');
        padp('layer_anchor','anchor');
        p.parse(varargin{:});
        obj.setParams(p.Results);
        
        
    end
    function run(obj,isfast)
        if nargin < 2
            isfast = false;
        end
        obj.elements = {};
        obj.addNBs(isfast);
        obj.addTether();
        obj.addExpoBoxes();
        obj.addCoupler();
        obj.addHolder();
        if obj.P_undercutMsk.isGen
            obj.addUndercutMsk();
        end
        obj.addZigZag();
        % move such that the buffer lower edge (or coupler endpoint) is at y = 0
        dy = - (obj.l_NB_unscaled/2 + obj.d_nb2tether + obj.d_tether2edge);
        obj.translate(0, dy);
        obj.vs_EL = obj.vs_EL + [0,0;dy,dy];
        obj.vs_ER = obj.vs_ER + [0,0;dy,dy];
    end
    function P = getDefaultUndercutMskParams(obj)
        P.d = 3;      % mask edge to NB edge distance
        P.w = 62.5;     % width of region above coupler
        P.h = 20;       % height of region above coupler
        P.isGen = false;
        P.layer = 'undercutMsk';
    end
    function P = getDefaultTetherParams(obj)
        P.w1 = 0.8;         % waveguide width, in um
        P.w2 = 2.4;         % waveguide width in connection with tether
        P.L = 10;           % width changing taper region length
        P.w_tether = 0.5;   % tether width
%         P.layer = 'M1_Neg';
    end
    function P = getDefaultCouplerParams(obj)
        P.L_taper = 10;     % taper length from NB width to WG width
        P.w = 0.8;          % waveguide width
        P.N_remove = 0;
        P.N_mirror = 18;
%         P.holeLayer = 'M1_Pos_holes';
%         P.slabLayer = 'M1_Neg';
        P.N_coupler = 8;
        P.isfullcav = 0;
        P.l_extra = 10;     % extra coupler length into the buffer region
    end
    function P = getDefaultZigzagParams(obj)
        P.w_tether = 0.1;
        P.w = 0.4;      % zigzag width
        P.g = 0.3;      % zigzag gap
        P.g_metal = 0.1;% gap between metal
        P.L = 5;        % zigzag arm length
        P.L2 = 20;      % zigzag total length
        P.d_NB2ZZ = 5;  % distance of NB to zigzag
        P.layer = 'zigzag';
        P.layer_metalSB = 'metal';
        P.layer_metalBB = 'metal_BB';
    end
end
methods (Access = private)
    function addCoupler(obj)
        obj.P_coupler.P_defect = obj.P_mirror;
        obj.P_coupler.P_mirror = obj.P_mirror;
        obj.P_coupler.N_trans = 0;
        obj.P_coupler.holeLayer = obj.layer_holes;
        obj.P_coupler.slabLayer = obj.layer_neg;
        l_remove = obj.P_coupler.N_remove * obj.P_coupler.P_mirror.a;  % removed length
        l_tot = obj.P_coupler.N_mirror * obj.P_coupler.P_mirror.a;
        h1 = obj.l_NB_unscaled + obj.d_nb2tether + obj.d_tether2edge;
        % waveguide with width P_coupler.w startes after the
        % P_coupler.L_taper transition region and ends at the edge between
        % exposure box 1 and box 2 (box 2 = buffer box for dicing)
        x_tetherEnd = obj.l_NB_unscaled/2 + obj.d_nb2tether + obj.P_tether.L;
        l_tetherEnd2Edge = obj.d_tether2edge - obj.P_tether.L;
        
        % generate coupler Polygon        
        % assuming no displacement or rotation of nanobeam
        dy = (obj.l_NB_unscaled/2 - obj.P_mirror.a/2)*obj.scale_NB;
        % right-bottom corner of NB before 90 deg rotation:
        x0 = (l_tot - obj.P_mirror.a/2) * obj.scale_NB - dy; y0 = obj.P_mirror.w/2 * obj.scale_NB;
        if isfield(obj.P_mirror,'isFBH')
            y0 = (obj.P_mirror.w - 2 * obj.P_mirror.amp)/2 * obj.scale_NB;
        end
        Ls = [obj.P_coupler.L_taper, x_tetherEnd - (x0 + obj.P_coupler.L_taper) ,...
            obj.P_coupler.l_extra + l_tetherEnd2Edge];
        Ws = [obj.P_coupler.w, obj.P_coupler.w, obj.w_end];
        obj.P_coupler.Ls = [];
        obj.P_coupler.Ws = Ws;
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
        couplerSlab = Polygon(coordCell,obj.layer_neg); 
        couplerSlab.rotate(90); 
        obj.addelement(couplerSlab); 
            
        coupler = nanobeamcoupler(rmfield(obj.P_coupler,{'L_taper','w','l_extra'}));
        coupler.zoom(obj.scale_NB);
        coupler.rotate(90); 
        % dunno why need an extra a/2; extra a/2 moved to l_tot
        coupler.translate(0, - dy );
        obj.addelement(coupler); obj.ind_coupler = length(obj.elements);
        
    end
    function nb = genNB(obj, isfast)
        nb = nanobeam('P_defect',obj.P_defect,'P_mirror',obj.P_mirror,...
            'N_trans',obj.N_trans,'N_mirror',obj.N_mirror,'isfullcav',1,...
            'holeLayer',obj.layer_holes,'slabLayer', obj.layer_neg,'isuptodate',false);
        if ~isempty(obj.maker)
            nb.maker = obj.maker;
        end
        nb.gennb(isfast);
        nb.zoom(obj.scale_NB);
        nb.rotate(90,0,0);
    end
    function addNBs(obj, isfast)
        nb1 = obj.genNB(isfast);
        nb1.translate(-obj.d_nb2refl2,0);%%%
        obj.addelement(nb1);  obj.ind_NB1 = length(obj.elements);
        nb2 = obj.genNB(isfast);
        nb2.translate(obj.d_nb2refl,0);
        obj.addelement(nb2);  obj.ind_NB2 = length(obj.elements);
                
        obj.l_NB_unscaled = obj.elements{obj.ind_NB1}.L;        
    end
    function addTether(obj)
        buff = 0.1;
        expandStructure(obj.P_tether);
        
        vTR = [w1/2; L];     vMR = [w2/2; 0];     vBR = [w1/2; -L];
        vBL = [-w1/2; -L];   vML = [-w2/2; 0];    vTL = [-w1/2; L];

        wg_taper = Polygon({vTR, vMR, vBR, vBL, vML, vTL});        
        w_box = (obj.d_nb2refl+obj.d_nb2refl2)  + 2 * obj.d_nb2edge;%%%
        tether = Rect(0,0,w_box + buff, w_tether,'base','center');
        % add extra triangle for holder strength
        % further add rect region beyond triangle
        dx = w_box/4; xMax = dx + w_box/2;
        dy = w_box/2/sqrt(3);
        w_rect = 5;
        tri1 = Polygon({[dx;0],[xMax;dy],[xMax + w_rect;dy],...
            [xMax + w_rect;-w_tether/2],[xMax;-w_tether/2]});
        tri2 = Polygon({[-dx;0],[-xMax;dy],[-xMax - w_rect;dy],...
            [-xMax - w_rect;-dy],[-xMax;-dy]});
        r_anchor = 15/2;
        c1 = Circ([xMax;r_anchor-w_tether],r_anchor);
        c1.layer = obj.layer_anchor;
        c2 = Circ([-xMax; 0],r_anchor);
        c2.layer = obj.layer_anchor;
        
        gtether = gpack.Group(0,obj.l_NB_unscaled/2 + obj.d_nb2tether,...
            {wg_taper, tether, tri1, tri2, c1, c2});
        gtether.layer = obj.layer_neg;
        obj.addelement(gtether);
        
        % copy to make HSQ cure layer
        g_cure = gtether.copy();
        g_cure.layer = obj.layer_cure;
        obj.addelement(g_cure);
    end
    function addZigZag(obj)
        % add two zigzags to the right nanobeam
        
        buff = 0.2;
        l_NB = obj.l_NB_unscaled * obj.scale_NB;
        x_zz = obj.d_nb2refl + obj.P_zigzag.d_NB2ZZ;
        y_zz = l_NB/2 + obj.P_mirror.w/2;   % holder width is the same as NB, see obj.addHolder
        zz1 = obj.genZigZag();
        zz1.translate(x_zz, y_zz);
        obj.vs_EL = obj.vs_EL + [x_zz;y_zz];
        obj.vs_ER = obj.vs_ER + [x_zz;y_zz];
        zz2 = obj.genZigZag();
        zz2.translate(x_zz, -y_zz);
        obj.vs_EL = obj.vs_EL + [0,x_zz;0,-y_zz];
        obj.vs_ER = obj.vs_ER + [0,x_zz;0,-y_zz];
        obj.addelement(zz1);
        obj.addelement(zz2);
        % add tether between zigzags and nanobeam
        w_tether = obj.P_zigzag.w_tether;
        l_tether = obj.P_mirror.w;
        r1 = Rect(x_zz - l_tether/2, y_zz - buff - w_tether/2, l_tether, w_tether, 'base','center');
        r1.layer = obj.layer_neg;
        r2 = Rect(x_zz - l_tether/2, -y_zz - buff + w_tether/2, l_tether, w_tether, 'base','center');
        r2.layer = obj.layer_neg;
        obj.addelement(r1);
        obj.addelement(r2);
    end
    function dev = genZigZag(obj)
        % This function includes parameters of the zig zag bender
        % param is a struct with field 'w', 'g' and 'L2'
        % 
        % enlarged tether, WTJ, 20180808
        % 
        % Modified from genZZBender_longArm_181002.m, WTJ, 20181016
        %
        param = obj.P_zigzag;
        w = param.w;
        g = param.g;

        zzp.w = w;
        zzp.L = param.L;
        zzp.L2 = param.L2;
        zzp.g = g;
        zzp.g_metal = param.g_metal;
        zzp.n = round(zzp.L2/(zzp.w+zzp.g)/2)+0.5;   % estimate the number of periods
        zzp.L2 = zzp.n*(zzp.w+zzp.g)*2 - zzp.w;      % calculate actual L2

        zz_LN = zigzag('w',zzp.w, 'L', zzp.L, 'g1', zzp.g, 'g2', zzp.g, 'n', zzp.n);

        zzp.w_metal = (zzp.w-zzp.g_metal)/2;
        zzp.g2 = zzp.g + zzp.w_metal * 2 + zzp.g_metal * 2;
        zz_M1 = zigzag('w',zzp.w_metal, 'L', zzp.L, 'g1', zzp.g2, 'g2', zzp.g, 'n', zzp.n);
        zz_M1.translate([-zzp.w/2 + zzp.w_metal/2;0]);
        zz_M2 = zigzag('w',zzp.w_metal, 'L', zzp.L, 'g1', zzp.g, 'g2', zzp.g2, 'n', zzp.n);
        zz_M2.translate([zzp.w/2 - zzp.w_metal/2;0]);

        zz_LN.layer = param.layer;
        zz_M1.layer = param.layer_metalSB;
        zz_M2.layer = param.layer_metalSB;
        % add tether    
        y0 = zz_LN.port1(2,1);
        x1 = zz_LN.port1(1,1);
        x2 = zz_LN.port1(1,2);
        dy = 6;
        x3 = x2 + dy;
        x4 = x1 - dy;
        y1 = y0 - dy;
        y2 = y1 - dy;
        p_tether = Polygon({[x2;y0],[x3;y1],[x3;y2],[x4;y2],[x4;y1],[x1;y0]});
        p_tether.layer = obj.P_zigzag.layer;
        r_rotCenter = [(x1+x2)/2;y0];
        p_tether.translate(-zzp.w/2, zzp.w/2);
        p_tether.rotate(r_rotCenter, 90);
        % add extra disk for tether/anchor
        r = 15/2;
        c_anchor = Circ((x1+x2)/2,y0 - 6-3, r, obj.layer_anchor);
        c_anchor.rotate(r_rotCenter, 90);
        g_sharps = obj.genSharpOnArc([(x1+x2)/2;y0 - 6-3], r, -225, 0);
        g_sharps.rotate(r_rotCenter, 90);
        % add BB region and release window
        w_BB = zzp.L + zzp.w * 2 + 1;
        h_BB = zzp.L2 + zzp.w * 2 + 1;
        r_BB = Rect(0,0,w_BB, h_BB, 'base','center');
        r_BB.layer = 'cure';
        r_release = Rect(0,0,w_BB + 3 + 3, h_BB, 'base','center');
%         r_release.layer = 'M2_undercutMsk';
        r_release.layer = obj.P_undercutMsk.layer;
        
        % add electrodes
        port_e1 = zz_M1.port1;      % end coords of the left electrode
        port_e2 = zz_M2.port1;      % end coords of the right electrode
        % left electrode extension
        v1 = port_e1(:,1); v2 = port_e1(:,2); v0 = [v1(1); v1(2) - zzp.w];
        v3 = [v2(1); v2(2) - zzp.w/2];
        v_endEL = [v1(1)+2*dy - 2.5;v1(2)-dy]; 
        v4 = [v_endEL(1);v_endEL(2) + zzp.w/2]; v5 = [v_endEL(1);v_endEL(2) - zzp.w/2];
        p_EL = Polygon({v0,v1,v2,v3,v4,v5});
        w_pad = 3;
%       r_EL = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
        c_EL = Circ(v_endEL(1), v_endEL(2), w_pad);
        g_EL = gpack.Group(0,0,{p_EL,c_EL});
        g_EL.layer = param.layer_metalSB;
%         r_EL_BB = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
%         r_EL_BB.layer = param.layer_metalBB;
        c_EL_BB = Circ(v_endEL(1), v_endEL(2), w_pad);
        c_EL_BB.layer = param.layer_metalBB;
        % right electrode extension
        v1 = port_e2(:,1); v2 = port_e2(:,2); v0 = [v1(1); v1(2) + zzp.w];
        v_endER = [v1(1)+2*dy;v1(2)+dy]; 
        v3 = [v_endER(1); v_endER(2) - zzp.w/2];
        v4 = [v_endER(1); v_endER(2) + zzp.w/2];
        p_ER = Polygon({v0,v1,v2,v3,v4});
        w_pad = 3;
%         r_ER = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
        c_ER = Circ(v_endER(1), v_endER(2), w_pad);
        g_ER = gpack.Group(0,0,{p_ER,c_ER});
        g_ER.layer = param.layer_metalSB;
%         r_ER_BB = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
%         r_ER_BB.layer = param.layer_metalBB;
        c_ER_BB = Circ(v_endER(1), v_endER(2), w_pad);
        c_ER_BB.layer = param.layer_metalBB;
        
        if strcmpi(param.layer_metalBB,'trash')
            g_EL.layer = 'trash';    % no big square
            g_ER.layer = 'trash';    % no big square
        end
        
        % make zigzag topleft at (0,0)
%         disp_zigzag = [param.L/2;-param.L2/2];
        disp_zigzag = [param.L/2;-zz_LN.port2(2,1)];
        dev = gpack.Group(disp_zigzag(1),disp_zigzag(2),...
            {zz_LN, zz_M1, zz_M2,p_tether,r_BB,r_release,c_anchor,g_sharps,...
            g_EL, c_EL_BB, g_ER, c_ER_BB});
        obj.vs_EL = [obj.vs_EL (v_endEL + disp_zigzag)];
        obj.vs_ER = [obj.vs_ER (v_endER + disp_zigzag)];
    
    end
    function g_sharps = genSharpOnArc(obj, pos, r, theta1, theta2)
        % add sharpe feature on arc with center pos and radius r, starting
        % at theta1 and end at theta2, in degree
        
        % parameter for single rhombus fixed:
        prd = 0.6;
        w_rhom = 0.2;
        h_rhom = w_rhom / 2 / tand(15/2);
        x = w_rhom/2; y = h_rhom;
        % function that generates single rhombus
        gsr = @()(Polygon({[x;0],[0;y],[-x;0],[0;-y]}));
        dtheta = prd/r/pi*180;
        thetas = [theta1:dtheta:theta2];
        xs = pos(1) + r * cosd(thetas);
        ys = pos(2) + r * sind(thetas);
        g_sharps = gpack.Group(0,0,{});
        for ii = 1:length(thetas)
            p_rhom = gsr();
            p_rhom.rotate(thetas(ii)+90);    % the rhombus is perpto radial direction
            p_rhom.translate(xs(ii), ys(ii));
            g_sharps.addelement(p_rhom);
        end
        g_sharps.layer = obj.layer_holes;
    end
    function addExpoBoxes(obj)
        % add exposure region
%         l_NB = l_NB_unscaled * obj.scale_NB;
        w1 = (obj.d_nb2refl+obj.d_nb2refl2) + 2 * obj.d_nb2edge; %%%
        h1 = obj.l_NB_unscaled + obj.d_nb2tether + obj.d_tether2edge;
        w2 = obj.w_buffer; 
        h2 = obj.h_buffer;
        b1 = Rect(0, h1/2-obj.l_NB_unscaled/2, w1,h1, 'base','center');
        b1.layer = obj.layer_pos;
        obj.addelement(b1);
        if w2 > 0 && h2 > 0
            b2 = Rect(0,h1 - obj.l_NB_unscaled/2 + h2/2, w2,h2,'base','center');
            b2.layer = obj.layer_pos;
            obj.addelement(b2);
        end
        % writefield box
        w_wf = 62.5;
        y_couplerEnd = obj.l_NB_unscaled/2 + obj.d_nb2tether + obj.d_tether2edge;
        rect_wf = Rect(0, y_couplerEnd - w_wf/2, w_wf, w_wf,...
            'base','center');
        rect_wf.layer = obj.layer_wf;
        obj.addelement(rect_wf);
    end
    function addUndercutMsk(obj)
        % add XeF2 undercut mask
        yMin = - obj.l_NB_unscaled * obj.scale_NB /2;
        yMax = obj.l_NB_unscaled/2 + obj.d_nb2tether + obj.d_tether2edge;
        % box width for NB region
        w = (obj.d_nb2refl+obj.d_nb2refl2) + obj.P_mirror.w * obj.scale_NB + obj.P_undercutMsk.d * 2;%%%
        % box width in coupler region
        w_coupler = obj.w_end + obj.P_undercutMsk.d * 2;
        h = obj.l_NB_unscaled;
        h_coupler = yMax - yMin - h;
        r1 = Rect(0, yMin + h/2, w, h, 'base','center');
        r2 = Rect(0, yMin + h + h_coupler/2, w_coupler, h_coupler, 'base','center');
        gMask = gpack.Group(0,0,{r1,r2});
        gMask.layer = obj.P_undercutMsk.layer;
        % copy for cure region
        g_cure = gMask.copy();
        g_cure.layer = obj.layer_cure;
        obj.addelement(g_cure);
        
        % add rect for dicing buffer region
        if obj.P_undercutMsk.h > 0
            r3 = Rect(0, yMax + obj.P_undercutMsk.h/2, obj.P_undercutMsk.w,...
                obj.P_undercutMsk.h, 'base','center');
            gMask.addelement(r3);
        end
        obj.addelement(gMask);
    end
    function addHolder(obj)
        % add holder for the NBs
        buff = 0.1;
        % first two top holders
        l1 = obj.d_nb2edge + obj.P_mirror.w + buff;
        w1 = obj.P_mirror.w;    % use same width as the NB
        dx = obj.d_nb2refl2 + obj.d_nb2edge/2 + buff/2;%%%
        dy = obj.l_NB_unscaled * obj.scale_NB/2 + w1/2 - buff;
        r2 = Rect(-dx, dy, l1, w1, 'base','center');
        % for right NB:
        l1_r = obj.P_zigzag.d_NB2ZZ;
        dx_r = obj.d_nb2refl + l1_r/2 - obj.P_mirror.w/2;
        r1 = Rect(dx_r, dy, l1_r, w1, 'base','center');        
        r1_bot = Rect(dx_r, -dy, l1_r, w1, 'base','center');
        % add extra triangles for holder strength, hard coded!
        % further add rect region beyond triangle
        dx_tri = 8; xMax = dx + dx_tri;
        yMin = dy - dx_tri/sqrt(3); yMax = dy + dx_tri/sqrt(3);
        w_rect = 5;
%         tri1 = Polygon({[dx;dy],[xMax;yMax],[xMax + w_rect;yMax],...
%             [xMax + w_rect;yMin],[xMax;yMin]});
        tri2 = Polygon({[-dx;dy],[-xMax;yMax],[-xMax - w_rect;yMax],...
            [-xMax - w_rect;yMin],[-xMax;yMin]});
        r_anchor = 15/2;
        c2 = Circ([-xMax;dy],r_anchor);
        c2.layer = obj.layer_anchor;
%         topHolder = gpack.Group(0,0,{r1,r2,tri1,tri2});
        topHolder = gpack.Group(0,0,{r1,r2,tri2, c2, r1_bot});
        topHolder.layer = obj.layer_neg;
        
        obj.addelement(topHolder);
        % copy for HSQ curing
        g_cure = topHolder.copy();
        g_cure.layer = obj.layer_cure;
        obj.addelement(g_cure);
        
        % one bottom holders, compensate for NB scaling
        w2 = (obj.d_nb2refl+obj.d_nb2refl2)  + 2 * obj.d_nb2edge;%%%
        h2 = 60 * 0.1;  % should be sufficient even for 60 um length NB with 90% scaling
        dy = - obj.l_NB_unscaled * obj.scale_NB /2 - h2/2 + buff;
%         rBot = Rect(0, dy, w2, h2, 'base','center');
%         rBot.layer = obj.layer_neg;
        % use polygon instead of rect to reduce area, WTJ, 180615
        yTop = dy + h2/2; yBot = dy - h2/2;
        x1 = obj.d_nb2refl2 + obj.P_mirror.w / 2 * obj.scale_NB;
        x2 = x1 + h2;
        h_rect = 5;
        x_max = obj.P_mirror.w/2;    % max right x coordinate for bottom tether, at gap middle
        pBot = Polygon({[x_max;yTop],[x_max;yBot - h_rect],...
            [-x2;yBot - h_rect],[-x2;yBot],[-x1;yTop]});
        pBot.layer = obj.layer_neg;
        rBot_pos = Rect(0, dy + h2, w2, h2, 'base','center');
        rBot_pos.layer = obj.layer_pos;
        
        c2 = Circ([x_max - r_anchor + 3; yBot - h_rect/2],r_anchor);
        c2.layer = obj.layer_anchor;
        
        obj.addelement(pBot);
        obj.addelement(rBot_pos);
        obj.addelement(c2);
        % copy for HSQ curing
        p_cure = pBot.copy();
        p_cure.layer = obj.layer_cure;
        obj.addelement(p_cure);
        
    end
    function setParams(obj, P)
        % set parameters
        % obj need to have properties as fields in P
        flds = fieldnames(P);
        for ii = 1:length(flds)
            fld = flds{ii};
            obj.(fld) = P.(fld);
        end
    end
end
end