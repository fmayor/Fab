function [dev,v_endEL,v_endER,v_endELN] = genSmallAnchorElectrodes(param)
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
    coeff_x3=param.coeff_x3;
    coeff_y3=param.coeff_y3; 
    h_elec_link=param.h_elec_link;
    g=param.g;
    dev = gpack.Group(0,0,{});

    layer_metalSB='metal';
    layer_metalBB='metal_BB';
    
    % middle electrode extension
    v1 = port_L(:,2);v5=port_L(:,1);
    v_endEL = [v1(1)+coeff_x*1;v1(2)-coeff_y*1];%dy
    v2=v1+[-0.1;-0.6];v3=v_endEL;v4=v3+[-0.8;0];
    p_EL = Polygon({v1,v2,v3,v4,v5});

%     r_EL = Rect(v_endEL(1), v_endEL(2), w_pad, w_pad,'base','center');
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
    v1 = port_R(:,1);
    v_endER = [v1(1)+coeff_x2*1;v1(2)-coeff_y2*1];%dy
    v6=v1;v7=v6+[0;-0.7];v9=v7+[-0.2;-0.3];v8=v_endER;
    v11=v6+[zz_g_metal;0];
    v12=[v11(1)+0.4;v8(2)];
    p_ER = Polygon({v6,v7,v9,v8,v12,v11});
    
    %r_ER = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    c_ER = Circ(v_endER(1), v_endER(2),w_pad);
    g_ER = gpack.Group(0,0,{p_ER,c_ER});
    g_ER.layer = layer_metalSB;
    dev.addelement(g_ER);
    %r_ER_BB = Rect(v_endER(1), v_endER(2), w_pad, w_pad,'base','center');
    %r_ER_BB.layer = layer_metalBB;
    c_ER_BB = Circ(v_endER(1), v_endER(2), w_pad);
    c_ER_BB.layer = layer_metalBB;
    dev.addelement(c_ER_BB);
    
    % left electrode extension
    v13 = port_R(:,1)+[zz_g_metal+(zz_w-zz_g_metal)/2;h_elec_link+zz_w+g/2];
    v18 = port_R(:,2)+[zz_g_metal+(zz_w-zz_g_metal)/2;h_elec_link+zz_w+g/2];
    v_endELN = [v13(1)+coeff_x3*1;v13(2)-coeff_y3*1];%dy
    v14=v13+[0.1;-1.5];
    v15=v_endELN;v16=v15+[1;0];v17=v14+[0.4;0];
    p_ELN = Polygon({v13,v14,v15,v16,v17,v18});

    %r_ELN = Rect(v_endELN(1), v_endELN(2), w_pad, w_pad,'base','center');
    c_ELN = Circ(v_endELN(1), v_endELN(2), w_pad);
    g_ELN = gpack.Group(0,0,{p_ELN,c_ELN});
    g_ELN.layer = layer_metalSB;
    dev.addelement(g_ELN);
%     r_EL_BBN = Rect(v_endELN(1), v_endELN(2), w_pad, w_pad,'base','center');
%     r_EL_BBN.layer = layer_metalBB;
    c_EL_BBN = Circ(v_endELN(1), v_endELN(2), w_pad);
    c_EL_BBN.layer = layer_metalBB;
    dev.addelement(c_EL_BBN);
    
    if param.is_mirrored
        dev.mirror([0;v1(2)],[1;0]);
        dev.mirror([v1(1);0],[0;1]);
        dev.translate([-(zz_w-zz_g_metal)/2;0]);
    end
end