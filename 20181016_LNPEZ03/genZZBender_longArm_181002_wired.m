function [dev,v_endEL,v_endER] = genZZBender_longArm_181002_wired(param)
% This function includes parameters of the zig zag bender
% param is a struct with field 'w', 'g' and 'L2'
% 
% enlarged tether, WTJ, 20180808
% 
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

    zz_LN.layer = 'M1_Neg';
    zz_M1.layer = 'metal';
    zz_M2.layer = 'metal';
    layer_metalSB = param.layer_metalSB;
    layer_metalBB = param.layer_metalBB;
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
    p_tether.layer = 'M1_Neg';
    r_rotCenter = [(x1+x2)/2;y0];
    p_tether.translate(0, zzp.w/2/sqrt(2));
    p_tether.rotate(r_rotCenter,45);
    % add extra disk for tether/anchor
    r = 15/2;
    c_anchor = Circ((x1+x2)/2,y0 - 6-3, r,'M1_LD');
    c_anchor.rotate(r_rotCenter,45);
    %% add sharps on anchor to help metal climb wall
    g_sharps = genSharpOnArc(param,[(x1+x2)/2;y0 - 6-3], r, -180, 0);
    g_sharps.rotate(r_rotCenter, 45);
    %% add metal electrodes on anchor
    
    port_e1 = zz_M1.port1;      % end coords of the left electrode
    port_e2 = zz_M2.port1;      % end coords of the right electrode
    % left electrode extension
    v1 = port_e1(:,1); v2 = port_e1(:,2); v0 = [v1(1); v1(2) - zzp.w];
    v3 = [v2(1); v2(2) - zzp.w/2];
    v_endEL = [v1(1)+0.75*dy;v1(2)-2.25*dy];
    v4 = [v_endEL(1);v_endEL(2) + zzp.w/2]; v5 = [v_endEL(1);v_endEL(2) - zzp.w/2];
    p_EL = Polygon({v0,v1,v2,v3,v4,v5});
    w_pad = 4;
    r_EL = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
    g_EL = gpack.Group(0,0,{p_EL,r_EL});
    g_EL.layer = layer_metalSB;
    r_EL_BB = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
    r_EL_BB.layer = layer_metalBB;
    % right electrode extension
    v1 = port_e2(:,1); v2 = port_e2(:,2); v0 = [v1(1); v1(2) + zzp.w/2];
    v_endER = [v1(1)+2.25*dy;v1(2)-0.75*dy];
    v3 = [v_endER(1); v_endER(2) - zzp.w/2];
    v4 = [v_endER(1); v_endER(2) + zzp.w/2];
    p_ER = Polygon({v0,v1,v2,v3,v4});
    w_pad = 4;
    r_ER = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    g_ER = gpack.Group(0,0,{p_ER,r_ER});
    g_ER.layer = layer_metalSB;
    r_ER_BB = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    r_ER_BB.layer = layer_metalBB;
        
    %% add BB region and release window
    w_BB = zzp.L + zzp.w * 2 + 1;
    h_BB = zzp.L2 + zzp.w * 2 + 3.5;
    r_BB = Rect(0,0,w_BB, h_BB, 'base','center');
    r_BB.layer = 'cure';
    r_release = Rect(0,0,w_BB + 3 + 3, h_BB, 'base','center');
    r_release.layer = 'M2_undercutMsk';    
    
%     %% Big electrodes
%     
%     w_pad = 200;
%     w_noLOR = 190;
%     w_wire = 5;
%     w_total = 62.5*4;
%     y_pad1 = -300;
% 
% 
% 
% 
%     d_wire2pad = (w_total-w_pad)/4;
%     y_pad2 = y_pad1 - w_pad - d_wire2pad*2;
%     w_pad3 = 260;
%     y_pad3 = y_pad2 - w_pad /2 - w_pad3/2- d_wire2pad*2;
%     y_line1 = y_pad2 - w_pad/2 - d_wire2pad;
%     y_line2 = y_pad3;
% 
% 
%     % wires from device to pads
%     % left electrodes connect to pad1 and pad2, right electrodes connect to line1
%     if param.flip==false
%         p1 = genBondPad(param,-v_endEL(1),y_pad1,w_pad);
%         p2 = genBondPad(param,-v_endEL(1),y_pad2,w_pad);
%         p3 = genBondPad(param,-v_endEL(1),y_pad3,w_pad3);
%         rect_line1 = Rect(-v_endEL(1), y_line1, w_total, w_wire,'base','center');
%         rect_line2 = Rect(-v_endEL(1),y_line2, w_total, w_wire,'base','center');
% 
%         % wire from ZZ2EL to pad1
%         v1 = v_endEL; v2 = [-v_endEL(1);y_pad1]; x1 = v1(1);
%         wire_ZZ2EL = Wire({v1,v2},w_wire);
%         % wire from ZZ2ER to line1
%         x2 = max(x1+ d_wire2pad,w_pad/2 + d_wire2pad);
%         v1 = v_endER; v2 = [50-v_endEL(1); v1(2)];
%         v3 = [50-v_endEL(1); y_pad1/2];
%         v4 = [x2-v_endEL(1); y_pad1/2];
%         v5 = [x2-v_endEL(1); (y_pad2 - w_pad /2 - d_wire2pad-w_wire/2)]; %v3 = [x2; (y_pad1+y_pad2)/2];
%         %v4 = [-w_pad/2 - d_wire2pad, (y_pad1+y_pad2)/2];
%         %v5 = [-w_pad/2 - d_wire2pad, y_line1];
%         wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);%wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);
%         g_wires = gpack.Group(0,0,{rect_line1, rect_line2, wire_ZZ2EL, wire_ZZ2ER});%wire_ZZ1EL, wire_ZZ1ER
%         g_wires.layer = layer_metalBB;
%         g_BB = gpack.Group(0,0,{p1,p2,p3, g_wires});
%     else
%         % wire from ZZ2EL to pad1
%         v1 = v_endEL; v2 = [v_endEL(1);y_pad1/2]; x1 = v1(1);
%         x2 = max(x1+ d_wire2pad,w_pad/2 + d_wire2pad);
%         v3 = [x2-v_endEL(1);y_pad1/2];
%         v4 = [x2-v_endEL(1);(y_pad1+y_pad2)/2];
%         v5 = [-x2-v_endEL(1)-w_wire/2;(y_pad1+y_pad2)/2];
%         wire_ZZ2EL = Wire({v1,v2,v3,v4,v5},w_wire);
%         % wire from ZZ2ER to line1
%         x2 = max(x1+ d_wire2pad,w_pad/2 + d_wire2pad);
%         x3=x2+7;
%         v1 = v_endER; v2 = [50-v_endEL(1); v1(2)];
%         v3 = [50-v_endEL(1);y_pad1/4];
%         v4 = [x3-v_endEL(1);y_pad1/4];
%         v5 = [x3-v_endEL(1); y_pad2]; %v3 = [x2; (y_pad1+y_pad2)/2];
%         v6 = [-v_endEL(1); y_pad2];
%         %v5 = [-w_pad/2 - d_wire2pad, y_line1];
%         wire_ZZ2ER = Wire({v1,v2,v3,v4,v5,v6},w_wire);%wire_ZZ2ER = Wire({v1,v2,v3,v4,v5},w_wire);
%         
%         g_wires = gpack.Group(0,0,{wire_ZZ2EL, wire_ZZ2ER});%wire_ZZ1EL, wire_ZZ1ER
%         g_wires.layer = layer_metalBB;
%         g_BB = gpack.Group(0,0,{g_wires});
%     end
% %     % wire from ZZ1EL to pad2
% %     x3 = x2 + d_wire2pad;
% %     v1 = vs_EL(:,1); v2 = [x3; v1(2)]; v3 = [x3; y_pad2]; v4 = [0;y_pad2];
% %     wire_ZZ1EL = Wire({v1,v2,v3,v4},w_wire);
% %     % wire from ZZ1ER to line1
% %     x4 = x3 + d_wire2pad;
% %     v1 = vs_ER(:,1); v2 = [x4; v1(2)]; v3 = [x4; y_line1];
% %     wire_ZZ1ER = Wire({v1,v2,v3},w_wire);




    
    
    
    %%
    dev = gpack.Group(0,0,{zz_LN, zz_M1, zz_M2,p_tether,r_BB,r_release,...
                            c_anchor,g_EL, r_EL_BB, g_ER, r_ER_BB,g_sharps});%,g_BB
    
    
    
end

function g_sharps = genSharpOnArc(param,pos, r, theta1, theta2)
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
        g_sharps.layer = param.layer_sharp;
end
    
function dev = genBondPad(param,x,y, w_bondPad)

if nargin < 3
    w_bondPad = 190;
end
    w_bondPadWindow = w_bondPad - 10;
    r_bondPad = Rect(0, 0, w_bondPad, w_bondPad, 'base','center');
    r_bondPad.layer = param.layer_metalBB;
    r_bondPadWin = Rect(0, 0, w_bondPadWindow, w_bondPadWindow, 'base','center');
    r_bondPadWin.layer = param.layer_uct;
    dev = gpack.Group(x,y,{r_bondPad,r_bondPadWin});
end
