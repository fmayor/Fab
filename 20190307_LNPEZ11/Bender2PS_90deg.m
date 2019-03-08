classdef Bender2PS_90deg < gpack.Group
   % a single-arm bender attached to a 1D phononic shield (PS)
   % Bender and PS parallel to x direction, from left to right is anchor,
   % bender, 1DPS
   % the last PS UC will be centered at (0,0).
   % the whole structure is symmetric for y-reflection
   %
   % WTJ, 20190225
   properties
       P_1DPS           % 1D PS parameters. See obj.getDefault1DPS for fields
       r_anchor         % radius of anchor pad
       P_bender         % bender parameters
       layer_LD         % low dose layer name
       layer_HD         % high dose layer name
       layer_RM         % layer name for release mask
       layer_MT         % layer name for metal
       layer_MT_PL      % layer name for photolitho metal
       layer_cure
       w_RM             % width of release mask
       r_metal          % radius of overlapping region between M2 and M3 metal
       cfg_wiring       % wiring config, 'T' for mid and top, 'B' for mid and bottom
       w_tip            % width of a small tip for touching
       l_tip            % length of a small tip, including
       dx_tip           % extra displacement for the tip
   end
   
   methods
       function obj = Bender2PS_90deg()
           % no input param allowed
           %
            obj@gpack.Group(0,0,{});
            obj.P_1DPS = obj.getDefault1DPSParams();
            obj.P_bender = obj.getDefaultBenderParams();
            obj.layer_HD = 'M1_HD';       % layer name of high dose
            obj.layer_LD = 'M1_LD';       % layer name of low dose
            obj.layer_MT = 'M2';       % layer name of metal
            obj.layer_MT_PL = 'M3';       % layer name of metal
            obj.layer_RM = 'M4';       % layer name of release mask    
            obj.layer_cure = 'M1_cure';  
            obj.w_RM = 5;
            obj.r_anchor = 8;
            obj.r_metal = 3;
            obj.cfg_wiring = 'T';
            obj.w_tip = 0;
            obj.l_tip = 0;
            obj.dx_tip = 0;
       end
       
       function run(obj)
           obj.elements = {};
           obj.add1DPS();
           obj.addBender();
           obj.addAnchor();
           obj.addReleaseMask();
           obj.addWiring();
           obj.addTip();
           obj.addProbePads();
       end
       
       function add1DPS(obj)
           g_1DPS = gpack.Group(0,0, {});
           g_1DPS.layer = obj.layer_HD;
           
           a = obj.P_1DPS.a;
           N = obj.P_1DPS.N;
           l_1DPS = obj.get1DPSLength();
           x_1DPS = - (l_1DPS - a)/2;
           r_tether = Rect(x_1DPS, 0, l_1DPS, obj.P_1DPS.w2,'base','center');
           g_1DPS.addelement(r_tether);
           
           xs = [-a*(N-1):a:0];
           for ii = 1:length(xs)
               x = xs(ii);
               g_1DPS.addelement(...
                    Rect(x, 0, obj.P_1DPS.w1l, obj.P_1DPS.w1w, 'base','center'));
           end
           if obj.cfg_wiring == 'T'
               g_1DPS.rotate(90);
               g_1DPS.translate([0;-0.725]);
           else
               g_1DPS.rotate(-90);
               g_1DPS.translate([0;0.725]);
           end
           
           obj.addelement(g_1DPS);
       end
       
       function addBender(obj)
           x_BenderRightEnd = obj.P_1DPS.w1w/4;%-obj.get1DPSLength() + obj.P_1DPS.a/2 ;
           
           L = obj.P_bender.L;
           w = obj.P_bender.w;
           w_metal = obj.P_bender.w_metal;
           g_metal = obj.P_bender.g_metal;
           
           r_LN = Rect(0,0, L, w, 'base','center');
           r_LN.layer = obj.layer_HD;
           y_metal = g_metal/2 + w_metal/2;
           r_metalT = Rect(0, y_metal, L, w_metal, 'base','center');
           r_metalT.layer = obj.layer_MT;
           r_metalB = Rect(0, -y_metal, L, w_metal, 'base','center');
           r_metalB.layer = obj.layer_MT;
           g_bender = gpack.Group(x_BenderRightEnd - L/2,0,...
               {r_LN, r_metalB, r_metalT});
           obj.addelement(g_bender);
           
       end
       
       function addAnchor(obj)
           x_buff = 0.1;
           x_anchor =  obj.P_1DPS.w1w/4 -...
               obj.P_bender.L - obj.r_anchor + x_buff;
           c_anchor = Circ([x_anchor;0], obj.r_anchor);
           c_anchor.layer = obj.layer_LD;
           obj.addelement(c_anchor);
       end
       
       function addReleaseMask(obj)
           % left margin
           x_L =  obj.P_1DPS.w1w/4 -...
               obj.P_bender.L;
           % right margin
           x_R = obj.w_RM;
           r_RM = Rect((x_L+x_R)/2, 0, x_R - x_L, 2*obj.w_RM,...
               'base','center');
           r_RM.layer = obj.layer_RM;
           obj.addelement(r_RM);
           % use the same for curing/cure
           r_cure = Rect((x_L+x_R)/2, 0, x_R - x_L, 2*obj.w_RM,...
               'base','center');
           r_cure.layer = obj.layer_cure;
           obj.addelement(r_cure);
       end
       
       function addWiring(obj)           
           x_anchor =  obj.P_1DPS.w1w/4 -...
               obj.P_bender.L - obj.r_anchor;
           w_final = 2;
           
           % port center coordinates based on wiring config
           % first port is wired to the top electrode of the bender
           if obj.cfg_wiring == 'T'
               x1 = x_anchor;
               y1 = obj.r_anchor;
               x2 = x_anchor - obj.r_anchor;
               y2 = 0;
           else               
               x1 = x_anchor - obj.r_anchor;
               y1 = 0;
               x2 = x_anchor;
               y2 = -obj.r_anchor;
           end
           x_benderLeft = x_anchor + obj.r_anchor;
           w_bdr = obj.P_bender.w;
           w_metal = obj.P_bender.w_metal;
           % left coordinates of bender electrodes
           vt1 = [x_benderLeft; w_bdr/2];
           vt2 = [x_benderLeft; w_bdr/2 - w_metal];
           vb1 = [x_benderLeft; -w_bdr/2 + w_metal];
           vb2 = [x_benderLeft; -w_bdr/2];
                      
           v_p1_1 = [x1; y1 + w_final/2];
           v_p1_2 = [x1; y1 - w_final/2];
           v_p2_1 = [x2; y2 + w_final/2];
           v_p2_2 = [x2; y2 - w_final/2];
           % polygon connecting between bender electrodes and pads
           p1 = Polygon({v_p1_1, vt1, vt2, v_p1_2});
           p2 = Polygon({v_p2_1, vb1, vb2, v_p2_2});
           pad1 = obj.genMetalPad();
           pad1.translate(x1,y1);
           pad2 = obj.genMetalPad();
           pad2.translate(x2,y2);
           g_wiring = gpack.Group(0,0,{p1,p2, pad1,pad2});
           g_wiring.layer = obj.layer_MT;
           
           obj.addelement(g_wiring);           
       end
       
       function g_MP = genMetalPad(obj)
           % generate petal pad between EBeam metal and photolitho metal
           c1 = Circ([0;0], obj.r_metal);
           c1.layer = obj.layer_MT;
           c2 = Circ([0;0], obj.r_metal);
           c2.layer = obj.layer_MT_PL;
           g_MP = gpack.Group(0,0, {c1,c2});
       end
       
       function addTip(obj)
           if obj.w_tip < 0.01
               return;
           end
           r_tip = Rect(obj.dx_tip, 0, obj.w_tip, obj.l_tip, 'base','center');
           r_tip.layer = obj.layer_HD;
           if obj.cfg_wiring == 'T'
               r_tip.translate([0;-obj.get1DPSLength()-...
                   obj.P_bender.w/2-obj.l_tip/2]);
           else
               r_tip.translate([0;obj.get1DPSLength()+...
                   obj.P_bender.w/2+obj.l_tip/2]);
           end
           obj.addelement(r_tip);
       end
       
       function addProbePads(obj)
           x_anchor =  obj.P_1DPS.w1w/4 -...
               obj.P_bender.L - obj.r_anchor;
           xL = x_anchor - obj.r_anchor;
           
           % y-position of the shrinked pads
            dy_shrinked = 10;   % y-separation of shrinked pads
            w_shrink = 6;
            y_pads = 0;
           if obj.cfg_wiring == 'T'
               y_pads = dy_shrinked/2;
           else
               y_pads = - dy_shrinked/2;
           end
           g_pads = obj.genProbePads();
           g_pads.translate(xL, y_pads);
           obj.addelement(g_pads);
           
           % the probe pads near y = 0 should already be connected
           % now connecting pads near y = +- r_anchor
           xR = x_anchor;
           if obj.cfg_wiring == 'T'
               yL = dy_shrinked;
               yR = obj.r_anchor;
           else
               yL = -dy_shrinked;
               yR = -obj.r_anchor;
           end
           wL = w_shrink;
           wR = obj.r_metal * 2;
           p = Polygon({[xL; yL + wL/2], [xR; yR + wR/2],...
               [xR; yR - wR/2], [xL, yL - wL/2]});
           p.layer = obj.layer_MT_PL;
           obj.addelement(p);
       end
       
       
       %% 
       function P = getDefault1DPSParams(obj)
           P.w1w = 0.7;
           P.w1l = 0.65;
           P.w2 = 0.05;
           P.a = 1;
           P.layer = 'M1_HD';
           P.N = 5;
       end
       function P = getDefaultBenderParams(obj)
           P.L = 10;
           P.w = 0.45;
           P.w_metal = 0.15;
           P.g_metal = 0.15;
       end
       function l = get1DPSLength(obj)
           l = obj.P_1DPS.N * obj.P_1DPS.a;
       end
       
        function setParams(obj, P)
            fns = fieldnames(P);
            for ii = 1:length(fns)
                fn = fns{ii};
                obj.(fn) = P.(fn);
            end
        end
        
        function g_pads = genProbePads(obj)
            % the center between the lower two shrinked pads are at (0, 0) for
            % cfg_wiring = 'T' and the center between the top two shrinked pads at
            % (0,0) for cfg_wiring = 'B'.
            w = 60;
            l = 100;
            dy = 100;
            if obj.cfg_wiring == 'T'
                y0 = -dy/2;
            else
                y0 = dy/2 - 2*dy;
            end
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
   end
    
end