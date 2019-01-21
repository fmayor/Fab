function dev = genCantilever(param)
% This function includes parameters of the zig zag bender
% param is a struct with field 'w', 'g' and 'L2'
% 
% enlarged tether, WTJ, 20190120
% this function is basically genLongArm with labels
% 
w = param.w;
g = param.g;

    zzp.w = w;
    zzp.L = param.L;
    zzp.L2 = param.L2;
    zzp.g = g;
    zzp.g_metal = 0.15;
    zzp.n = round(zzp.L2/(zzp.w+zzp.g)/2)+0.5;   % estimate the number of periods
    zzp.L2 = zzp.n*(zzp.w+zzp.g)*2 - zzp.w;      % calculate actual L2

    zz_LN = zigzag('w',zzp.w, 'L', zzp.L, 'g1', zzp.g, 'g2', zzp.g, 'n', zzp.n);

    zzp.w_metal = (zzp.w-zzp.g_metal)/2;
    zzp.g2 = zzp.g + zzp.w_metal * 2 + zzp.g_metal * 2;
    zz_M1 = zigzag('w',zzp.w_metal, 'L', zzp.L, 'g1', zzp.g2, 'g2', zzp.g, 'n', zzp.n);
    zz_M1.translate([-zzp.w/2 + zzp.w_metal/2;0]);
    zz_M2 = zigzag('w',zzp.w_metal, 'L', zzp.L, 'g1', zzp.g, 'g2', zzp.g2, 'n', zzp.n);
    zz_M2.translate([zzp.w/2 - zzp.w_metal/2;0]);

    zz_LN.layer = 'M1_MD'; %M1_Neg (changed on 20181106)
    zz_M1.layer = 'metal';
    zz_M2.layer = 'metal';
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
    % add BB region and release window
    w_BB = zzp.L + zzp.w * 2 + 1;
    h_BB = zzp.L2 + zzp.w * 2 + 1;
    r_BB = Rect(0,0,w_BB, h_BB, 'base','center');
    r_BB.layer = 'cure';
    r_release = Rect(0,0,w_BB + 3 + 3, h_BB, 'base','center');
    r_release.layer = 'M2_undercutMsk';    
    

    dev = gpack.Group(0,0,{zz_LN, zz_M1, zz_M2,r_BB,r_release,c_anchor,p_tether});
    
    %add small electrodes
    v1=[x1,y0];
    wire1 = Wire(x1,y0,0.15);
    wire1.layer = 'metal';
    dev.addelement(wire1);
    
    %add label
    [ind_pos, ind_neg] = genStringPolygon(sprintf('%d',param.ind),...
            0.8);
    ind_pos.layer = 'M1_Neg';
    ind_neg.layer = 'M1_Pos_ind';
    ind = gpack.Group(0,0,{ind_pos, ind_neg});
    dx = -15;
    dy = -20;
    ind.rotate(-90);
    ind.translate(dx,dy);

    dev.addelement(ind);
    
    %add writefield
    wf = Rect(0,0,62.5, 62.5, 'base','center');
    wf.translate(0,-15);
    wf.layer='M1_Field';
    dev.addelement(wf);
    
end
