classdef TPNC < gpack.Group
    % tunable phononic crystal (PNC) on Y cut 
    % composed of few PNCs and two Bender2PSs
    % The center of the PNC defect will be at x = 0
    % the top most PNC will be at y = 0
    %
    % Total size ~ 400 by 600 um
    %
    % WTJ, 20190225
    %
    
    properties
        P_PNC
        P_Bender2PS_90deg
        
        N_PNC           % # pf PNCs
        L_Bender        % bender length
        dy_PNC          % PNC crystals y separation
        x_tip_L         % left bender tip position
        x_tip_R         % right bender tip position
        
        g_PS2PNCL        % gap between PS and PNC for left tip
        g_PS2PNCR       % gap between PS and PNC for right tip
        layer_HD        % layer name of high dose
        layer_LD        % layer name of low dose
        layer_MT        % layer name of metal
        layer_MT_PL        % layer name of photolitho metal
        layer_RM        % layer name of release mask        
        layer_Field
        layer_Field2
        layer_cure
        
        Ps_PNC_def      % PNC defect parameters
        inds_PNC         % index for PNC elements
        ind_BenderL
        ind_BenderR
        
        
    end
    
    methods
        function obj = TPNC()
            % no param allowed
            
            obj@gpack.Group(0,0,{});
            obj.P_PNC = obj.getDefaultPNCParams();
            obj.P_Bender2PS_90deg = obj.getDefaultBender2PSParams();
            obj.dy_PNC = 4;
            obj.N_PNC = 5;
            P_PNC_def.w1w_def = 1;
            P_PNC_def.w1l_def = 1;
            obj.Ps_PNC_def = repmat(P_PNC_def, 1, obj.N_PNC);
            obj.inds_PNC = [];
            obj.ind_BenderL = [];
            obj.ind_BenderR = [];
            obj.g_PS2PNCL = 0.25;
            obj.g_PS2PNCR = 0.25;
            obj.layer_HD = 'M1_HD';       % layer name of high dose
            obj.layer_LD = 'M1_LD';       % layer name of low dose
            obj.layer_MT = 'M2';       % layer name of metal
            obj.layer_MT_PL = 'M3';       % layer name of metal
            obj.layer_RM = 'M4';       % layer name of release mask    
            obj.layer_Field = 'M1_Field';
            obj.layer_cure = 'M1_cure';
            obj.layer_Field2 = 'M2_Field';
            
            obj.L_Bender = 10;
            obj.x_tip_L = -0.1; %-1.2
            obj.x_tip_R = 1.25; %2.2
            
            
        end
        function run(obj)
            obj.elements = {};
            obj.syncParams();
            obj.addPNCs();
            obj.addBender2PSs();
            % generate port at the lowest PNC
            dy_PNCPort = -(obj.N_PNC - 1)*obj.dy_PNC;
            for ii = 1:obj.N_PNC
                obj.elements{obj.inds_PNC(ii)}.dy_port = dy_PNCPort + (ii-1)*obj.dy_PNC;
                obj.elements{obj.inds_PNC(ii)}.run();
                obj.elements{obj.inds_PNC(ii)}.translate(0, -(ii-1)*obj.dy_PNC);
            end
            % remove one cell from 1DPS of the bender if the touchpoint is two cells away
%             N_rem = round(abs(obj.x_tip_L /( 2 * obj.P_PNC.a)));
%             obj.elements{obj.ind_BenderL}.P_1DPS.N = obj.P_PNC.N_mirror - N_rem;
%             N_rem = round(abs(obj.x_tip_R /( 2 * obj.P_PNC.a)));
%             obj.elements{obj.ind_BenderR}.P_1DPS.N = obj.P_PNC.N_mirror - N_rem;
            
            obj.elements{obj.ind_BenderL}.run();
            obj.elements{obj.ind_BenderR}.run();
            obj.elements{obj.ind_BenderR}.rotate(180);
            %dy = obj.P_PNC.w1w/2 + obj.g_PS2PNC +...
                %obj.P_Bender2PS_90deg.l_tip/2;
            w_defs=[obj.Ps_PNC_def.w1w_def];
            dyL = obj.g_PS2PNCL+obj.P_Bender2PS_90deg.l_tip+0.175+0.45/2+...
                obj.P_Bender2PS_90deg.P_1DPS.N*obj.P_Bender2PS_90deg.P_1DPS.a+w_defs(1)/2;
            dyR = obj.g_PS2PNCR+obj.P_Bender2PS_90deg.l_tip+0.175+0.45/2+...
                obj.P_Bender2PS_90deg.P_1DPS.N*obj.P_Bender2PS_90deg.P_1DPS.a+obj.P_PNC.w1w/2;
            
            obj.elements{obj.ind_BenderL}.translate(obj.x_tip_L, dyL);
            obj.elements{obj.ind_BenderR}.translate(obj.x_tip_R, dyR);
            
            obj.addWriteField();
            obj.addProbePads();
        end
        
        function syncParams(obj)
            % layer names
            layernames = {'layer_HD','layer_LD','layer_MT',...
                'layer_MT_PL', 'layer_RM', 'layer_cure'};
            for ii = 1:length(layernames)
                fn = layernames{ii};
                obj.P_PNC.(fn) = obj.(fn);
                obj.P_Bender2PS_90deg.(fn) = obj.(fn);
            end
            obj.P_PNC.l_arm = obj.L_Bender;
            obj.P_Bender2PS_90deg.P_1DPS = obj.P_PNC;
            obj.P_Bender2PS_90deg.P_1DPS.N = 4;%obj.P_PNC.N_mirror;
        end
        
        function setN_PNC(obj, N)
            obj.N_PNC = N;            
            P_PNC_def.w1w_def = 1;
            P_PNC_def.w1l_def = 1;
            obj.Ps_PNC_def = repmat(P_PNC_def, 1, obj.N_PNC);
        end
        
        function addPNCs(obj)
            obj.inds_PNC = ones(1, obj.N_PNC);
            for ii = 1:obj.N_PNC
                pnc = PNC;
                pnc.setParams(obj.P_PNC);
                % override defect params:
                pnc.setParams(obj.Ps_PNC_def(ii));
                obj.addelement(pnc);
                obj.inds_PNC(ii) = length(obj.elements);
            end            
        end
        
        function addBender2PSs(obj)
            bdr2psL = Bender2PS_90deg;
            bdr2psL.setParams(obj.P_Bender2PS_90deg);
            bdr2psL.P_bender.L = obj.L_Bender;
            obj.addelement(bdr2psL);
            obj.ind_BenderL = length(obj.elements);
            
            bdr2psR = Bender2PS_90deg;
            bdr2psR.setParams(obj.P_Bender2PS_90deg);
            bdr2psR.P_bender.L = obj.L_Bender;
            obj.addelement(bdr2psR);
            obj.ind_BenderR = length(obj.elements);
            bdr2psR.cfg_wiring = 'B';
        end
        
        function addProbePads(obj)
            % add probe pads for PNCs
            y_pads = -200;
            g_pads = obj.genProbePads();
            dy_PNCPort = -(obj.N_PNC - 1)*obj.dy_PNC;
            y_wireEnd = dy_PNCPort - 20;
            % length of taper
            l_taper = 20;
            
            % at this point, PNCs should be already run
            l_PNC = obj.elements{obj.inds_PNC(1)}.getTotalLength();
            x_pads = -l_PNC/2;
            
            g_pads.translate(x_pads, y_pads);
            obj.addelement(g_pads);
            w = obj.elements{obj.inds_PNC(1)}.r_metal * 2;
            dy_shrinked = 10;
            w_shrink = 6;
            
            % generate a path bending up from the probe pads
            w_pos = dy_shrinked + w_shrink;
            w_neg = dy_shrinked - w_shrink;
%             paths = {};
%             path{1} = [x_pads; y_pads];
            r_bend = l_PNC/2;
            arc_pos = Arc_zero([0; 0], w_pos, r_bend, 90);
            arc_pos.layer = obj.layer_MT_PL;
            arc_neg = Arc_zero([0; 0], w_neg, r_bend, 90);
            arc_neg.layer = [obj.layer_MT_PL '_Neg'];
            g_arc = gpack.Group(0,0,{arc_pos, arc_neg});
            g_arc.rotate(-90);
            g_arc.translate(x_pads, y_pads);
            obj.addelement(g_arc);
            
            wire_pos = Wire({[0; y_pads + r_bend], [0; y_wireEnd - l_taper],...
                [0; y_wireEnd]}, [w_pos, w_pos, l_PNC]);
            wire_neg = Wire({[0; y_pads + r_bend], [0; y_wireEnd - l_taper],...
                [0; y_wireEnd]}, [w_neg, w_neg, l_PNC/w_pos*w_neg]);
            wire_pos.layer = obj.layer_MT_PL;
            wire_neg.layer = [obj.layer_MT_PL '_Neg'];
            
            obj.addelement(wire_pos);
            obj.addelement(wire_neg);
            
            % final wire from wireEnd to PNC ports
            w_start = (l_PNC - l_PNC/w_pos*w_neg)/2;    % width of positive wire end
            w_stop = obj.elements{obj.inds_PNC(1)}.r_metal * 2;
            x_start = l_PNC/w_pos*w_neg/2 + w_start/2;
            x_stop = l_PNC/2;
            p1 = Polygon({[x_start-w_start/2;y_wireEnd],[x_stop-w_stop/2; dy_PNCPort],...
                [x_stop+w_stop/2; dy_PNCPort], [x_start+w_start/2;y_wireEnd]});
            p2 = Polygon({[-x_start-w_start/2;y_wireEnd],[-x_stop-w_stop/2; dy_PNCPort],...
                [-x_stop+w_stop/2; dy_PNCPort], [-x_start+w_start/2;y_wireEnd]});
            p1.layer = obj.layer_MT_PL;
            p2.layer = obj.layer_MT_PL;
            obj.addelement(p1);
            obj.addelement(p2);
            
        end
        
        
        function g_pads = genProbePads(obj)
            % copied from Bender2PS
            w = 60;
            l = 100;
            dy = 100;
            y0 = dy/2 - 2*dy;
            ys = [y0:dy:(y0+2*dy)];
            g_pads = gpack.Group(0,0,{});
            for ii = 1:length(ys)
                y = ys(ii);
                rect_metal = Rect(0, y, l, w, 'base','center');
                rect_metal.layer = obj.layer_MT_PL;
                r_RM = Rect(0, y, l - 3*2, w - 3*2, 'base','center');
                r_RM.layer = obj.layer_RM;
                g_pads.addelement(r_RM);
                g_pads.addelement(rect_metal);                
            end
            
            % generate shrink frp 60 um to 6 um
            dy_shrinked = 10;
            w_shrink = 6;
            xL = l/2;       % right end of pads, left end of shrink region
            l_shrink = 40;
            xR = xL + l_shrink;
            y1 = dy/2; y2 = -dy/2;
            y1s = dy_shrinked/2; 
            y2s = -dy_shrinked/2;
            p1 = Polygon({[xL; y1 + w/2], [xR; y1s + w_shrink/2],...
                [xR; y1s - w_shrink/2], [xL; y1-w/2]});
            p2 = Polygon({[xL; y2 + w/2], [xR; y2s + w_shrink/2],...
                [xR; y2s - w_shrink/2], [xL; y2-w/2]});
            p1.layer = obj.layer_MT_PL;
            p2.layer = obj.layer_MT_PL;
            g_pads.addelement(p1);
            g_pads.addelement(p2);
            g_pads.translate(-xR, 0);
        end
        
        function addWriteField(obj)
            y = -obj.N_PNC * (obj.dy_PNC-1)/2;
            w = 62.5;
            rect = Rect(0, y, w, w, 'base','center');
            rect.layer = obj.layer_Field;
            rect2 = Rect(0, y, w, w, 'base','center');
            rect2.layer = obj.layer_Field2;
            rect3 = Rect(w, y, w, w, 'base','center');
            rect3.layer = obj.layer_Field2;
            rect4 = Rect(-w, y, w, w, 'base','center');
            rect4.layer = obj.layer_Field2;
            obj.addelement(gpack.Group(0,0,{rect, rect2,rect3,rect4}));
        end
        %%
        function P = getDefaultPNCParams(obj)
            P.N_mirror = 7;
            P.a = 1;
            P.w1w = 0.7;
            P.w1l = 0.65;
            P.w2 = 0.05;
            P.w1w_def = 1;
            P.w1l_def = 1;
            
        end
        function P = getDefaultBender2PSParams(obj)
            P.w_RM = 4;
            P.r_anchor = 8;
            P.r_metal = 3;
            P.w_tip = 0.3; %0.05
            P.l_tip = 0.4; %1.2
            P.dx_tip = 0;
        end
    end
    
    
end