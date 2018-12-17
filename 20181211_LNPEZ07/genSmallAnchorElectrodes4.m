function [dev,v_endEL,v_endER] = genSmallAnchorElectrodes4(param)
    dy=param.dy;
    w_pad=param.w_pad;
    port_L=param.port_L;
    port_R=param.port_R;
    zz_w=param.zz_w;
    
    dev = gpack.Group(0,0,{});

    layer_metalSB='metal';
    layer_metalBB='metal_BB';
    
    % left electrode extension
    v1 = port_L(:,1); v2 = port_L(:,2); v0 = [v1(1); v1(2) - 1];
    v3 = [v2(1); v2(2) - zz_w/2];
    v_endEL = [v1(1)+0.75*dy;v1(2)-2.25*dy];
    v4 = [v_endEL(1);v_endEL(2) + 0.5]; v5 = [v_endEL(1);v_endEL(2) - 0.5];
    p_EL = Polygon({v0,v1,v2,v3,v4,v5});
    %r_EL = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
    c_EL = Circ(v_endEL(1), v_endEL(2), w_pad);
    g_EL = gpack.Group(0,0,{p_EL,c_EL});
    g_EL.layer = layer_metalSB;
    dev.addelement(g_EL);
%     r_EL_BB = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
%     r_EL_BB.layer = layer_metalBB;
    c_EL_BB = Circ(v_endEL(1), v_endEL(2), w_pad);
    c_EL_BB.layer = layer_metalBB;

    dev.addelement(c_EL_BB);
    % right electrode extension
    v1 = port_R(:,1); v2 = port_R(:,2); v0 = [v1(1); v1(2) + zz_w/2];
    v_endER = [v1(1)+2.25*dy;v1(2)-0.75*dy];
    v3 = [v_endER(1); v_endER(2) - 0.5];
    v4 = [v_endER(1); v_endER(2) + 0.5];
    p_ER = Polygon({v0,v1,v2,v3,v4});
    w_pad = 4;
%     r_ER = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    c_ER = Circ(v_endER(1), v_endER(2), w_pad);
    g_ER = gpack.Group(0,0,{p_ER,c_ER});
    g_ER.layer = layer_metalSB;
    dev.addelement(g_ER);
    %r_ER_BB = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    %r_ER_BB.layer = layer_metalBB;
    c_ER_BB = Circ(v_endER(1), v_endER(2), w_pad);
    c_ER_BB.layer = layer_metalBB;
    dev.addelement(c_ER_BB);
end