function [dev,v_endEL,v_endER] = genSmallAnchorElectrodes3(param)
    %dy=param.dy;
    w_pad=param.w_pad;
    port_L=param.port_L;
    port_R=param.port_R;
    zz_w=param.zz_w;
    zz_g_metal = param.zz_g_metal;
    coeff_x=param.coeff_x;
    coeff_y=param.coeff_y;
    coeff_x2=param.coeff_x2;
    coeff_y2=param.coeff_y2;   
    dev = gpack.Group(0,0,{});

    layer_metalSB='metal';
    layer_metalBB='metal_BB';
    
    % left electrode extension
    v1 = port_L(:,1); v2 = port_L(:,2); v0 = [v1(1); v1(2) - zz_w];
    v3 = [v2(1); v2(2) - zz_w/2];
    v_endEL = [v1(1)+coeff_x*1;v1(2)-coeff_y*1];%dy
    v4 = [v_endEL(1);v_endEL(2) + zz_w/2]; v5 = [v_endEL(1);v_endEL(2) - zz_w/2];
%     v6= v1+[zz_g_metal/2;-0.3];v7=v6+[0;0.3];v8=[v_endEL(1);v7(2)];v9=[v_endEL(1);v_endEL(2)];
%     v10=v9+[zz_g_metal;0];v11=v9+[zz_g_metal;zz_g_metal];v12=v7+[zz_g_metal;0];
%     v13=v6+[zz_g_metal;0];
    v5 = v1+[0;+0.3]; v6=v5+[0;-0.4];v7=[v_endEL(1);v6(2)];v8=[v7(1);v_endEL(2)];
    v9=v8+[zz_g_metal;0];v10 = [v9(1);v6(2)-zz_g_metal];v11=[v2(1);v10(2)];
    v12=[v11(1);v5(2)];
    
    %v8=v1+[0;0.3];v9=v1;v10=v_endEL;v11=v_endEL+[0;zz_g_metal];v12=[v2(1);v11(2)];v13=v2+[0;0.3];
    %p_EL = Polygon({v0,v1,v2,v3,v4,v5});
    %p_EL = Wire({v1+[zz_g_metal/2;0.1],v_endEL},zz_g_metal);
    %p_EL = Wire({v6,v7,v_endEL},zz_g_metal);
    p_EL = Polygon({v5,v6,v7,v8,v9,v10,v11,v12});
    r_EL = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
    g_EL = gpack.Group(0,0,{p_EL,r_EL});
    dev.addelement(g_EL);
    g_EL.layer = layer_metalSB;
    r_EL_BB = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
    r_EL_BB.layer = layer_metalBB;
    dev.addelement(r_EL_BB);
    % right electrode extension
    v1 = port_R(:,1); v2 = port_R(:,2); v0 = [v1(1); v1(2) + zz_w/2];
    v_endER = [v1(1)+coeff_x2*1;v1(2)-coeff_y2*1];%dy
    v3 = [v_endER(1); v_endER(2) - zz_w/2];
    v4 = [v_endER(1); v_endER(2) + zz_w/2];
    v5 = v1+[0;+0.3]; v6=v5+[0;-1];v7=[v_endER(1);v6(2)];v8=[v7(1);v_endER(2)];
    v9=v8+[zz_g_metal;0];v10 = [v9(1);v6(2)-zz_g_metal];v11=[v2(1);v10(2)];
    v12=[v11(1);v5(2)];
    %p_ER = Polygon({v0,v1,v2,v3,v4});
    %p_ER = Wire({v1+[zz_g_metal/2;0.1],v_endER},zz_g_metal);
    p_ER = Polygon({v5,v6,v7,v8,v9,v10,v11,v12});
    r_ER = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    g_ER = gpack.Group(0,0,{p_ER,r_ER});
    g_ER.layer = layer_metalSB;
    dev.addelement(g_ER);
    r_ER_BB = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    r_ER_BB.layer = layer_metalBB;
    dev.addelement(r_ER_BB);
    if param.is_mirrored
        dev.mirror([0;v1(2)],[1;0]);
    end
end