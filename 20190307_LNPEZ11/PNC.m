classdef PNC < gpack.Group
    % class for a single phononic crystal cavity
    % The object is centered at (0, 0), with PNC parallel to x direction.
    % No input params for constructor
    % change parameters after constructing the object
    % use obj.run to generate actual components after setting the
    % parameters
    %
    % WTJ
    % 20190225
    properties
        N_mirror        % number of mirrors
        a               % pitch of the mirrors
        w1w             % mirror cell width ( along y direction)
        w1l             % mirror cell length (along the direction of the PNC)
        w2              % tether width between mirror cells
        w1w_def          % defect cell width
        w1l_def          % defect cell length
        w_RM            % release mask width
        d_metal         % metal electrode separation
        r_metal         % radius of overlapping region between M2 and M3 metal
        r_anchor        % radius of the anchor
        layer_HD        % layer name of high dose
        layer_LD        % layer name of low dose
        layer_MT        % layer name of metal
        layer_MT_PL        % layer name of photolitho metal
        layer_RM        % layer name of release mask        
        layer_cure      % curing layer
        w_arm           % extra arm for compatibility with bender
        l_arm           % extra arm for compatibility with bender
        dy_port         % y direction shift of the ports
    end
    
    methods
        function obj = PNC
            % constructor for PNC
            % No input parameter allowed
            
            obj@gpack.Group(0,0,{});
            P = obj.getDefaultParams();
            obj.setParams(P);
            
        end
        
        function run(obj)
            obj.elements = {};
            obj.addMirrors();
            obj.addDefect();
            obj.addAnchors();
            obj.addArms();
            obj.addMetal();
            obj.addReleaseMask();
        end
        
        function addMirrors(obj)
            l_def = obj.a - obj.w1l + obj.w1l_def;
            xs_mir = 0:obj.a:((obj.N_mirror-1) * obj.a);
            xs_mir = l_def/2 + obj.a/2 + xs_mir;
            xs_mir = [-flip(xs_mir), xs_mir];           % centers of mirror cells
            g_mir = gpack.Group(0,0,{});
            for ii = 1:length(xs_mir)
                x = xs_mir(ii);
                g_mir.addelement(...
                    Rect(x, 0, obj.w1l, obj.w1w, 'base','center'));
            end
            % add tether as one Rect
            l_tether = obj.getCavityLength();
            g_mir.addelement(Rect(0,0, l_tether, obj.w2, 'base','center'));
            g_mir.layer = obj.layer_HD;
            obj.addelement(g_mir);
        end
        
        
        function addDefect(obj)
            r_def = Rect(0,0, obj.w1l_def, obj.w1w_def,'base','center');
            r_def.layer = obj.layer_HD;
            obj.addelement(r_def);
        end
        
        
        function addArms(obj)
            x_arm = obj.getCavityLength()/2 + obj.l_arm/2;
            r_arm_R = Rect(x_arm, 0, obj.l_arm, obj.w_arm,'base','center');
            r_arm_L = Rect(-x_arm, 0, obj.l_arm, obj.w_arm,'base','center');
            g_LN = gpack.Group(0,0,{r_arm_R, r_arm_L});
            g_LN.layer = obj.layer_HD;
            obj.addelement(g_LN);
            
            % add metal wire on arms
            g_armWire = g_LN.copy();
            g_armWire.layer = obj.layer_MT;
            obj.addelement(g_armWire);
        end
        
        
        function addAnchors(obj)
            x_buff = 0.1;
            x_anchor = obj.getCavityLength()/2 + obj.l_arm + obj.r_anchor - x_buff;
            r_R = Circ([x_anchor;0], obj.r_anchor);
            r_L = Circ([-x_anchor;0], obj.r_anchor);
            r_R.layer = obj.layer_LD;
            r_L.layer = obj.layer_LD;
            obj.addelement(r_R);
            obj.addelement(r_L);
            
        end
        
        function addMetal(obj)
            l_def = obj.a - obj.w1l + obj.w1l_def;
            xs_mir = 0:obj.a:((obj.N_mirror-1) * obj.a);
            xs_mir = l_def/2 + obj.a/2 + xs_mir;
            xs_mir = [-flip(xs_mir), xs_mir];           % centers of mirror cells
            g_metal = gpack.Group(0,0,{});
            for ii = 1:length(xs_mir)
                x = xs_mir(ii);
                g_metal.addelement(...
                    Rect(x, 0, obj.w1l, obj.w1w, 'base','center'));
            end
            % add tether as one Rect
            l_tether = (obj.getCavityLength() - obj.w1l_def)/2;
            x_tether = obj.getCavityLength()/2 - l_tether/2;
            g_metal.addelement(Rect(x_tether,0, l_tether, obj.w2, 'base','center'));
            g_metal.addelement(Rect(-x_tether,0, l_tether, obj.w2, 'base','center'));
            g_metal.layer = obj.layer_MT;
            % add metal on defect
            w_metal = obj.w1w_def;
            l_metal = (obj.w1l_def - obj.d_metal)/2;
            if l_metal < 0.1
                warning('Metal on defect too small! Using 100 nm.');
                l_metal = 0.1;
                % preventing l_metal too small
            end
            x_metal_def = obj.w1l_def/2 - l_metal/2;
            g_metal.addelement(Rect(x_metal_def,0, l_metal, w_metal, 'base','center'));
            g_metal.addelement(Rect(-x_metal_def,0, l_metal, w_metal, 'base','center'));
            
            % add wire from arm end to anchor ends
            l_wire = obj.r_anchor * 2;
            w_wire = obj.w_arm;         % use same width as the arm
            w_wireEnd = 2;
            x_wire = obj.getTotalLength()/2 - l_wire/2;
            xL = x_wire - l_wire/2;
            xR = x_wire + l_wire/2;
            yR = obj.dy_port;
            
            
            p_wire_R = Polygon({[xL; w_wire/2], [xR; yR+w_wireEnd/2],...
                [xR; yR-w_wireEnd/2], [xL; -w_wire/2]});            
            p_wire_L = Polygon({[-xL; w_wire/2], [-xR; yR+w_wireEnd/2],...
                [-xR; yR-w_wireEnd/2], [-xL; -w_wire/2]});
            g_metal.addelement(p_wire_R);
            g_metal.addelement(p_wire_L);
            
            
            % add M2 to M3 connection
            x_circ = obj.getTotalLength()/2;
            c_metal_L = Circ([x_circ; obj.dy_port], obj.r_metal);
            c_metal_R = Circ([-x_circ; obj.dy_port], obj.r_metal);
            g_metal.addelement(c_metal_L);
            g_metal.addelement(c_metal_R);
            
            g_metal_PL = gpack.Group(0,0,{c_metal_L, c_metal_R});   % photolitho metal
            g_metal_PL.layer = obj.layer_MT_PL;
            
            obj.addelement(g_metal);
            obj.addelement(g_metal_PL);
        end
        
        
        function addReleaseMask(obj)
            l_RM = obj.getCavityLength() + 2*obj.l_arm;
            r_RM = Rect(0,0, l_RM, obj.w_RM,...
                'base','center');
            r_RM.layer = obj.layer_RM;
            obj.addelement(r_RM);
            % use the same for curing/cure
            r_cure = Rect(0,0, l_RM, obj.w_RM,...
                'base','center');
            r_cure.layer = obj.layer_cure;
            obj.addelement(r_cure);
        end
        
        
        %% auxillary functions
        
        function P = getDefaultParams(obj)
            P.N_mirror = 5;
            P.a = 1;
            P.w1w = 0.7;
            P.w1l = 0.65;
            P.w2 = 0.05;
            P.w1w_def = 1;
            P.w1l_def = 1;
            P.r_anchor = 8;
            P.w_RM = 4+4;
            P.d_metal = 0.5;
            P.r_metal = 3;
            P.layer_HD = 'M1_HD';       % layer name of high dose
            P.layer_LD = 'M1_LD';       % layer name of low dose
            P.layer_MT = 'M2';       % layer name of metal
            P.layer_MT_PL = 'M3';       % layer name of metal
            P.layer_RM = 'M4';       % layer name of release mask     
            P.layer_cure = 'M1_cure';       % layer name of release mask     
            P.w_arm = 1.5;
            P.l_arm = 1;              % short by default
            P.dy_port = 0;
        end
        
        function ltot = getCavityLength(obj)
            % get defect + mirror length
            l_def = obj.a - obj.w1l + obj.w1l_def;          % length of the defect
            ltot = obj.a * obj.N_mirror * 2 + l_def;
        end
        function ltot = getTotalLength(obj)
            % get total device length
            ltot = obj.getCavityLength() + obj.l_arm * 2 + obj.r_anchor * 4;
        end
        
        
        function setParams(obj, P)
            fns = fieldnames(P);
            for ii = 1:length(fns)
                fn = fns{ii};
                obj.(fn) = P.(fn);
            end
        end
    end
end