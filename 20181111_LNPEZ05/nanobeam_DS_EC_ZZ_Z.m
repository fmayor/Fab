classdef nanobeam_DS_EC_ZZ_Z < gpack.Group
% nanobeam with edge coupler for lense fiber
% lense fiber -> coupler -> reflector, two nanobeams side-coupled to the
% reflector
% Vertically oriented
% WTJ, 20180523
% 
% Add secondary mask for XeF2, change first mask to adapt to HSQ
% WTJ, 20180615
% 
% Modified for strain-release holder
% WTJ, 20180712
% 
% Add electrodes
% WTJ, 20180719
% 
% Add 1D phonon shield
% WTJ, 20180808
% 
% Add sharpe corners for increasing electrode connectivity
% WTJ, 20180829

properties
    w_end       % coupler end width
    w_buffer    % buffer region width
    h_buffer    % buffer region depth
    P_tether    % tether parameters
    P_coupler   % coupler parameters
    P_undercutMsk   % undercut mask parameters
    P_electrode % electrode parameters
    P_1DPS      % 1D phonon shield parameters
    P_sharp     % sharpe shape (triangles) at the boundary of the anchor pads
                % for better electrode connectivity
    P_mirror
    P_defect
    N_trans
    N_mirror
    maker       % nanobeam unitcell maker
    d_nb2refl    % center distance between nanobeam and reflector
    d_nb2edge   % distance from nanobeam to undercut box edge
    d_tether2edge   % vertical distance from tether to buffer region edge
    scale_NB    % scaling factor for nanobeams, also scales the mirrors of the reflector
    ind_NB1     % index of the nanobeam object
    ind_NB2
    ind_coupler
    ind_gtether
    ind_expoBox
    l_NB_unscaled
    layer_pos   % layer name of area to expose (no resist, to etch)
    layer_holes
    layer_neg   % layer name of area not to expose
    layer_wf    % write field layer
    layer_cure  % layer for HSQ curing
    layer_anchor  % layer name for large area anchor
    layer_medium % layer name for medium dose zigzags
    h_link
    w_link
    w %zizag arm width
    L %zigzag arm length
    L2 %zigzag length in direction transversal to arm
    g_metal %zigzag gap between electrodes
    g %zigzag gap between arms
end
methods
    function obj = nanobeam_DS_EC_ZZ_Z(varargin)
        obj@gpack.Group(0,0,{});
        P_tether = obj.getDefaultTetherParams();
        P_coupler = obj.getDefaultCouplerParams();
        P_undercutMsk = obj.getDefaultUndercutMskParams();
        P_electrode = obj.getDefaultElectrodeParams();
        P_1DPS = obj.getDefault1DPSParams();
        P_sharp = obj.getDefaultPSharp();
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
        padp('d_nb2edge',10,@isnumeric);
%         padp('d_nb2tether',5,@isnumeric);
        padp('d_tether2edge',7,@isnumeric);
        padp('scale_NB',1,@isnumeric);
        padp('P_tether',P_tether,@isstruct);
        padp('P_coupler',P_coupler,@isstruct);
        padp('P_undercutMsk',P_undercutMsk,@isstruct);
        padp('P_mirror',P_mirror,@isstruct);
        padp('P_defect',P_defect,@isstruct);
        padp('P_electrode',P_electrode,@isstruct);
        padp('P_1DPS',P_1DPS,@isstruct);
        padp('P_sharp',P_sharp,@isstruct);
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
        obj.x = 0; obj.y = 0;
        obj.addNBs(isfast);
        obj.addTether_Holder();
        obj.addExpoBoxes();
        obj.addCoupler();
        if obj.P_undercutMsk.isGen
            obj.addUndercutMsk();
        end
        if obj.P_electrode.isGen
            obj.addElectrodes();
        end
        % move such that the buffer lower edge (or coupler endpoint) is at y = 0
        dy = - (obj.l_NB_unscaled/2 + obj.d_tether2edge);
        obj.translate(0, dy);
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
        P.w_tetherArm = 30; % tether length
    end
    function P = getDefaultCouplerParams(obj)
        P.L_taper = 10;     % taper length from NB width to WG width
        P.w = 0.8;          % waveguide width
        P.N_remove = 0;
        P.N_mirror = 18;
        P.N_coupler = 8;
        P.isfullcav = 0;
        P.l_extra = 10;     % extra coupler length into the buffer region
    end
    function P = getDefaultElectrodeParams(obj)
        % default parameters for 
        P.dx = 15;      % x-aligned electrode separation (along the NB)
        P.dy = 0.1;     % y-aligned electrode separation (tranverse the NB)
        P.d_metal2defect = 7.5;      % distance from electrode to defect cell 
                                    % y-aligned electrodes
        P.layer.metal1 = 'metal1';
        P.layer.metal2 = 'metal2';
        P.layer.metal1_BB = 'metal1_BB';   % bigbeam layer
        P.layer.metal2_BB = 'metal2_BB';
        P.isGen = false;
        P.config = 'B';             % 'L', 'R' or 'B', left only, right only,
                                    % or both NBs have electrodes
    end
    function P = getDefault1DPSParams(obj)
        % default parameters for 1D phonon shield
        % x direction is direction parallel to NB (global x in COMSOL)
        % PS is connected along y direction
        P.w1x = 0.6;
        P.w1y = 0.65;
        P.w2 = 0.05;
        P.a = 1;
        P.isGen = false;
    end
    function P = getDefaultPSharp(obj)
        % get default P_sharp
        P.isGen = false;
        P.theta = 30;   % 30 deg angle
        P.period = 2;
        P.layer = '';
        P.dutycyc = 0.5;  % duty cycle is the width of the triangle (rhombus) over period
        
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
        h1 = obj.l_NB_unscaled + obj.d_tether2edge;
        % waveguide with width P_coupler.w startes after the
        % P_coupler.L_taper transition region and ends at the edge between
        % exposure box 1 and box 2 (box 2 = buffer box for dicing)
        x_tetherEnd = obj.l_NB_unscaled/2 + obj.P_tether.L;
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
        nb1.translate(-obj.d_nb2refl,0);
        obj.addelement(nb1);  obj.ind_NB1 = length(obj.elements);
        nb2 = obj.genNB(isfast);
        nb2.translate(obj.d_nb2refl,0);
        obj.addelement(nb2);  obj.ind_NB2 = length(obj.elements);
                
        obj.l_NB_unscaled = obj.elements{obj.ind_NB1}.L;        
    end
    function addTether_Holder(obj)
        buff = 1;
        expandStructure(obj.P_tether);
        
        vTR = [w1/2; L];     vMR = [w2/2; 0];     vBR = [w1/2; -L];
        vBL = [-w1/2; -L];   vML = [-w2/2; 0];    vTL = [-w1/2; L];

        wg_taper = Polygon({vTR, vMR, vBR, vBL, vML, vTL});        
        w_box = obj.d_nb2refl * 2 + 2 * obj.d_nb2edge;
%         w_tetherArm = 30;   % mvoed to obj.P_tether as a parameter
        if obj.P_1DPS.isGen
            tether = obj.gen1DPS(w_tetherArm);
        else
            tether = Rect(0,0,w_tetherArm, w_tether,'base','center');
            %Make rectangle to be substracted to unlink outer zipper edge from tether so it can be
            %attached to zigzag
            w_rect_sub=1;
            h_rect_sub=2.55*obj.P_mirror.a*obj.scale_NB;
            x_pos=-obj.d_nb2refl-w_rect_sub/2;
            y_pos=-w_tether/2-h_rect_sub/2;%-dev.d_tether2edge-dev.P_tether.w_tether/2-h_rect/2;
            rect_TL = Rect(x_pos,y_pos,w_rect_sub,h_rect_sub, 'base','center');
            rect_TL.layer = obj.layer_holes; 
            %rect_TR = Rect(-x_pos,y_pos,w_rect_sub,h_rect_sub, 'base','center');
            %rect_TR.layer = obj.layer_holes; 
        end
        tether_cure = Rect(0,0,w_tetherArm, w_tether + 3 + 3,'base','center');
        tether_cure.layer = obj.layer_cure;
        % add extra triangle for holder strength
        % further add rect region beyond triangle
        w_pad = 10;
        dx = w_tetherArm/2 - buff; xMax = dx + w_pad/2 * sqrt(3);
        dy = w_pad/2;
        w_rect = 5;
        tri1 = Polygon({[dx;0],[xMax;dy],[xMax + w_rect;dy],...
            [xMax + w_rect;-dy],[xMax;-dy]});
        tri2 = Polygon({[-dx;0],[-xMax;dy],[-xMax - w_rect;dy],...
            [-xMax - w_rect;-dy],[-xMax;-dy]});
        % add circular plate for larger Xactix tolerance
        x_anchor = w_tetherArm/2 + w_pad/2* sqrt(3);
        r_anchor = 6;%15/2;
        c1 = Circ([x_anchor;0],r_anchor);
        c1.layer = obj.layer_anchor;
        c2 = Circ([-x_anchor;0],r_anchor);
        c2.layer = obj.layer_anchor;
        
        gtether = gpack.Group(0, (obj.l_NB_unscaled + 2*obj.P_mirror.a/3) * obj.scale_NB/2,...
            {wg_taper, tether,rect_TL, tri1, tri2, tether_cure, c1, c2});
        
        if obj.P_sharp.isGen
            g_sharp1 = obj.genSharpOnArc([x_anchor;0], r_anchor, -180+90, 180-90);
            g_sharp1.layer = obj.P_sharp.layer;
            g_sharp2 = obj.genSharpOnArc([-x_anchor;0], r_anchor, 90, 360-90);
            g_sharp2.layer = obj.P_sharp.layer;
            gtether.addelement(g_sharp1);
            gtether.addelement(g_sharp2);
        end
        gtether.layer = obj.layer_neg;
        obj.addelement(gtether);
        obj.ind_gtether = length(obj.elements);
        
        rect_BL = Rect(x_pos,y_pos+w_tether+h_rect_sub,w_rect_sub,h_rect_sub, 'base','center');
        rect_BL.layer = obj.layer_holes; 
        %rect_BR = Rect(-x_pos,y_pos+w_tether+h_rect_sub,w_rect_sub,h_rect_sub, 'base','center');
        %rect_BR.layer = obj.layer_holes; 
        gholder = gpack.Group(0, - (obj.l_NB_unscaled + 2*obj.P_mirror.a/3) * obj.scale_NB/2,...
            {tether,rect_BL, tri1, tri2, tether_cure, c1, c2});        
        if obj.P_sharp.isGen
            gholder.addelement(g_sharp1);
            gholder.addelement(g_sharp2);
        end
        gholder.layer = obj.layer_neg;
        obj.addelement(gholder);        
        
        
        % copy to make HSQ cure layer
        g_cure = gtether.copy();
        g_cure.layer = obj.layer_cure;
        obj.addelement(g_cure);
        g_cure = gholder.copy();
        g_cure.layer = obj.layer_cure;
        obj.addelement(g_cure);
        
        %make tether link
        x_pos_link = x_pos+obj.P_mirror.w/4-4*obj.w_link;
        y_pos_link = y_pos-h_rect_sub/2+obj.h_link/2;
        %rect_help = Rect(x_pos_link,y_pos_link-obj.h_link/2-0.3/2,0.64,0.3, 'base','center');
        link_TL = Rect(x_pos_link,y_pos_link,obj.w_link,obj.h_link, 'base','center');
        %link_TR = Rect(-x_pos_link,y_pos_link,obj.w_link,obj.h_link, 'base','center');
        link_BL = Rect(x_pos_link,-y_pos_link,obj.w_link,obj.h_link, 'base','center');
        %link_BR = Rect(-x_pos_link,-y_pos_link,obj.w_link,obj.h_link, 'base','center');
        
        
        %make zigzags
        x_pos_ZZ=x_pos_link-obj.L/2;
        y_pos_ZZ=y_pos_link+obj.h_link/2+obj.L2/2+obj.g/2+obj.w/2;
        zz_TL = genZZBender_longArm_181002(obj);
        zz_TL.translate([x_pos_ZZ;y_pos_ZZ]);
        %zz_TR = genZZBender_longArm_181002(obj);
        %zz_TR.mirror([0;0],[1;0]);
        %zz_TR.translate([-x_pos_ZZ;y_pos_ZZ]);
        gZZT = gpack.Group(0, (obj.l_NB_unscaled + 2*obj.P_mirror.a/3) * obj.scale_NB/2,...
            {link_TL,zz_TL});
        gZZT.layer = 'M1_MD';
        obj.addelement(gZZT);
        zz_BL = genZZBender_longArm_181002(obj);
        zz_BL.mirror([0;0],[1;0]);
        zz_BL.translate([x_pos_ZZ;-y_pos_ZZ]);
        %zz_BR = genZZBender_longArm_181002(obj);
        %zz_BR.translate([-x_pos_ZZ;-y_pos_ZZ]);
        gZZB = gpack.Group(0, - (obj.l_NB_unscaled + 2*obj.P_mirror.a/3) * obj.scale_NB/2,...
        {link_BL,zz_BL});
        gZZB.layer = 'M1_MD';
        obj.addelement(gZZB);
    end
    function g_rhm = genRhombusAry(obj, rs)
        % generate rhombus array on multiple straight lines
        % rs is an 2 by m matrix containing coordinates of m points
        [~, m] = size(rs);
        g_rhm = gpack.Group(0,0,{});
        for ii = 1:(m-1)
            g_rhm_single = obj.genRhombusAry_single(rs(:,ii), rs(:,ii+1));
            g_rhm.addelement(g_rhm_single);
        end
        if isempty(obj.P_sharp.layer)
            obj.P_sharp.layer = obj.layer_holes;
        end
        g_rhm.layer = obj.P_sharp.layer;
    end
    function g_rhm = genRhombusAry_single(obj, r_start, r_stop)
        % generate an array of rhombus, given the start and stop point
        % coordinates r_start and r_stop (i.e., two column vector)
        x1 = r_start(1); y1 = r_start(2);
        x2 = r_stop(1); y2 = r_stop(2);
        d = norm(r_start - r_stop);
        
        theta = atand((y2-y1)/(x2-x1));   % angle from x-axis
        prd = obj.P_sharp.period;
        w_rhom = prd * obj.P_sharp.dutycyc;
        h_rhom = w_rhom / 2 / tand(obj.P_sharp.theta/2);
        x = w_rhom/2; y = h_rhom;
        % function that generates single rhombus
        gsr = @()(Polygon({[x;0],[0;y],[-x;0],[0;-y]}));
        
        % generate array
        g_rhm = gpack.Group(0,0,{});
        n = floor(d/prd);
        xs = linspace(x1,x2,n);
        ys = linspace(y1,y2,n);
        for ii = 1:n
            rhm = gsr();
            rhm.rotate(theta);
            rhm.translate(xs(ii), ys(ii));
            g_rhm.addelement(rhm);
        end
    end    
    function g_sharps = genSharpOnArc(obj, pos, r, theta1, theta2)
        % add sharpe feature on arc with center pos and radius r, starting
        % at theta1 and end at theta2, in degree
        
        % parameter for single rhombus fixed:
        prd = obj.P_sharp.period;
        w_rhom = prd * obj.P_sharp.dutycyc;
        h_rhom = w_rhom / 2 / tand(obj.P_sharp.theta/2);
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
    function g_PS = gen1DPS(obj, w_arm)
        % generate 1D phonon shield
        w1_hori = obj.P_1DPS.w1y;
        w1_verti = obj.P_1DPS.w1x;
        w2 = obj.P_1DPS.w2;
        a = obj.P_1DPS.a;
        g_PS = gpack.Group(0,0,{});
        r = Rect(0,0,w_arm, w2,'base','center');
        % add rect for center part to hold NBs and reflectors
        w_tether = obj.P_tether.w_tether;
        w_center = obj.d_nb2refl * 2 + obj.P_mirror.w + obj.P_mirror.amp * 2;
        r_center = Rect(0,0,w_center, w_tether, 'base','center');
        g_PS.addelement(r);
        g_PS.addelement(r_center);
        % position of the bulk blocks
        xs = [-w_arm/2:a:w_arm/2];
        for ii = 1:length(xs)
            g_PS.addelement(Rect(xs(ii),0,w1_hori, w1_verti,...
                'base','center'));
        end
    end
    function addExpoBoxes(obj)
        % add exposure region
%         l_NB = l_NB_unscaled * obj.scale_NB;
        w1 = obj.d_nb2refl * 2 + 2 * obj.d_nb2edge;
        h1 = obj.l_NB_unscaled + obj.d_tether2edge;
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
        rect_wf = Rect(0, -obj.l_NB_unscaled * obj.scale_NB/2 - 7 + w_wf/2, w_wf, w_wf,...
            'base','center');
        rect_wf.layer = obj.layer_wf;
        obj.addelement(rect_wf);
        obj.ind_expoBox = length(obj.elements);
    end
    function addUndercutMsk(obj)
        % add XeF2 undercut mask
        yMin = - obj.l_NB_unscaled * obj.scale_NB /2;
        yMax = obj.l_NB_unscaled/2 + obj.d_tether2edge;
        % box width for NB region
        w = 2 * obj.d_nb2refl + obj.P_mirror.w * obj.scale_NB + obj.P_undercutMsk.d * 2;
        % box width in coupler region
        w_coupler = obj.w_end + obj.P_undercutMsk.d * 2;
        tol = 10;    % extra tolerance +-5 um
        h = obj.l_NB_unscaled * obj.scale_NB;
        h_coupler = yMax - yMin - h;
        r1 = Rect(0, yMin + h/2, w, h + tol, 'base','center');
        r2 = Rect(0, yMin + h + h_coupler/2, w_coupler, h_coupler, 'base','center');
        gMask = gpack.Group(0,0,{r1,r2});
        gMask.layer = obj.P_undercutMsk.layer;
        % copy for cure region
        
        w1_cure = 2 * obj.d_nb2refl + obj.P_mirror.w * obj.scale_NB + 4 * 2;
        w2_cure = 4 * 2;
        r1_cure = Rect(0, yMin + h/2, w1_cure, h, 'base','center');
        r2_cure = Rect(0, yMin + h + h_coupler/2, w2_cure, h_coupler, 'base','center');
        g_cure = gpack.Group(0,0,{r1_cure,r2_cure});
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
    function addElectrodes(obj)
        % add two types of electrodes
        % needs to boolean AND with the LN domain from first mask
        config = obj.P_electrode.config;
        [g_metal1, g_metal2] = obj.genLeftElectrodes();
        if config == 'L'        % left only
            obj.addelement(g_metal1);
            obj.addelement(g_metal2);
        elseif config == 'R'    % right only
            g_metal1.mirror([0;0],[0;1]);
            g_metal2.mirror([0;0],[0;1]);
            obj.addelement(g_metal1);
            obj.addelement(g_metal2);
        elseif config == 'B'    % both NBs
            [g_metal1_R, g_metal2_R] = obj.genLeftElectrodes();
            g_metal1_R.mirror([0;0],[0;1]);
            g_metal2_R.mirror([0;0],[0;1]);
            obj.addelement(g_metal1);
            obj.addelement(g_metal2);
            obj.addelement(g_metal1_R);
            obj.addelement(g_metal2_R);
            l_NB = (obj.l_NB_unscaled + 2*obj.P_mirror.a/3)* obj.scale_NB;
            % short the left and right electrodes near tether (upper)
            w_arm = 30;
            tol = 0.04;      % +- 20 nm alignment tolerance
            h_arm = obj.P_1DPS.w2 + tol;
            r_topArmMid = Rect( 0, l_NB/2,...
                w_arm, h_arm,'base','center');
            r_botArmMid = Rect( 0, -l_NB/2,...
                w_arm, h_arm,'base','center');
            if ~ isfield(obj.P_electrode.layer, 'metal1_abs')
                r_topArmMid.layer = obj.P_electrode.layer.metal1;
                r_botArmMid.layer = obj.P_electrode.layer.metal1;
            else
                r_topArmMid.layer = obj.P_electrode.layer.metal1_abs;
                r_botArmMid.layer = obj.P_electrode.layer.metal1_abs;
            end
            obj.addelement(r_topArmMid);
            obj.addelement(r_botArmMid);
        end
        
    end
    function [g_metal1, g_metal2] = genLeftElectrodes(obj)
        
        dx = obj.P_electrode.dx;
        dy = obj.P_electrode.dy;
        d_metal2defect = obj.P_electrode.d_metal2defect;
        l_NB = obj.l_NB_unscaled * obj.scale_NB;
        d_nb2refl = obj.d_nb2refl;
        
        % x-aligned electrodes
        w = 1.5;
        l = (l_NB - dx)/2 + 2;  % extra 2 um for safety
        r_top = Rect(- d_nb2refl, dx/2 + l/2, w, l, 'base','center');
        r_bot = Rect(- d_nb2refl, -dx/2-l/2, w, l, 'base','center');
        % extra area to connect to LN pads
        if isfield(obj.P_tether,'w_tetherArm')
            w_arm = obj.P_tether.w_tetherArm + 5;
        else 
            w_arm = 30;
        end
        h_arm = 5;
        r_topArm = Rect(- d_nb2refl - w_arm/2, l_NB/2,...
            w_arm, h_arm,'base','center');
        r_botArm = Rect(- w_arm/2, -l_NB/2,...
            w_arm, h_arm,'base','center');
        w_BB = 13;
        r_topBB = Rect(- w_arm, l_NB/2, w_BB, w_BB, 'base','center');
        r_topBB.layer = obj.P_electrode.layer.metal1_BB;
        r_botBB = Rect(- w_arm, -l_NB/2, w_BB, w_BB, 'base','center');
        r_botBB.layer = obj.P_electrode.layer.metal1_BB;
        g_metal1 = gpack.Group(0,0,{r_top, r_bot, r_topArm, r_botArm,...
            r_topBB, r_botBB});
        g_metal1.layer = obj.P_electrode.layer.metal1;
        
        % y-aligned electrodes
        w = 0.7;
        a_mir = obj.P_mirror.a * obj.scale_NB;
        l = (l_NB - 2*a_mir/3 - 2*d_metal2defect)/2;
        r_left = Rect(- d_nb2refl - dy/2 - w/2, -(l_NB-2*a_mir/3)/2 + l/2, w, l, 'base','center');
        r_right = Rect(- d_nb2refl + dy/2 + w/2, -(l_NB-2*a_mir/3)/2 + l/2, w, l, 'base','center');
        % bottom holder parameters (y position and width):
        y_holder = - (obj.l_NB_unscaled + 2*obj.P_mirror.a/3) * obj.scale_NB/2;
        w_tether = obj.P_tether.w_tether;
        r_right_extra = Rect(- d_nb2refl + dy/2 + w/2, y_holder, w, w_tether,...
            'base','center');
        % connecting fromNB to pad
        d_onHolder = w_tether/3;
        w_arm = 30;
        h_arm = 5;
        r_topHalf = Rect(- d_nb2refl - dy/2 - w_arm/2, y_holder + d_onHolder/2 + h_arm/2,...
            w_arm, h_arm, 'base','center');
        r_botHalf = Rect(- d_nb2refl + dy/2 + w - w_arm/2, y_holder - d_onHolder/2 - h_arm/2,...
            w_arm, h_arm, 'base','center');
        % bigbeam region
        w_BB = 5;
        r_topEnd = Rect(-25, y_holder + 5, w_BB, w_BB, 'base','center');
        r_topEnd.layer = obj.P_electrode.layer.metal2_BB;
        r_botEnd = Rect(-25, y_holder - 5, w_BB, w_BB, 'base','center');
        r_botEnd.layer = obj.P_electrode.layer.metal2_BB;
        
        g_metal2_SB = gpack.Group(0,0,{r_left, r_right,r_right_extra, r_topHalf,...
            r_botHalf});
        g_metal2_SB.layer = obj.P_electrode.layer.metal2;
        g_metal2 = gpack.Group(0,0,{g_metal2_SB,r_topEnd,r_botEnd});
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