function devAll = addBendedEdgeCoupler(dev, param)
% add a bended edge coupler to an existing nanobeam_DS_EC device
% Using hard coded parameters so far.
% WTJ, 20180919
global config;
R = 25;
w = dev.P_coupler.w;
w_end = dev.w_end;    % this is already after SEM correction
if isfield(param,'theta')
    theta = param.theta;
else
    theta = 45;
end

devAll = gpack.Group(0,0,{});


if theta == 0    
    % don't do anything, just add index
    dev.run(config.isfast);
%     ind = param.ind;
%     g_ind = genLNNB_DS_EC_ind(ind);
%     dev.addelement(g_ind);
    devAll.addelement(dev);
    return;
end

%% remove original end-coupler (after tether)
L_halfTether = dev.P_tether.L;
L_endCoupler = dev.d_tether2edge - L_halfTether;
dev.d_tether2edge = L_halfTether;
dev.w_end = dev.P_coupler.w;

% generate dev
dev.run(config.isfast);

% add index
% ind = param.ind;
% g_ind = genLNNB_DS_EC_ind(ind);
% dev.addelement(g_ind);
devAll.addelement(dev);
%% generate arc
arc = Arc_zero([0;0], w, R, theta);
arc.layer = dev.layer_neg;
arc_cure = Arc_zero([0;0], w + 6, R, theta);
arc_cure.layer = dev.layer_cure;
arc_undercut = Arc_zero([0;0], w + 8, R, theta);
arc_undercut.layer = dev.P_undercutMsk.layer;
devAll.addelement(arc);
devAll.addelement(arc_cure);
devAll.addelement(arc_undercut);

%% generate tether and end-coupler
% tether
theta_genExtraTether = 30;
if abs(theta) > theta_genExtraTether
    gtether = dev.elements{dev.ind_gtether};
    gtether_EC = gtether.copy();
    gtether_EC.y = 0;
    if dev.P_electrode.isGen
        % add electrodes for coupler holder
        d_elec2mid = 2.5;    % electrode edge to center
        l_elec = 18;
        w_elec = 5;
        x_elec = l_elec/2 + d_elec2mid;
        r_elec1 = Rect(x_elec,0,l_elec, w_elec,'base','center');
        r_elec1.layer = dev.P_electrode.layer.metal1;
        r_elec2 = Rect(-x_elec,0,l_elec, w_elec,'base','center');
        r_elec2.layer = dev.P_electrode.layer.metal1;
        r_elec_abs = Rect(0, 0, l_elec, 0.1, 'base','center');
        r_elec_abs.layer = dev.P_electrode.layer.metal1_abs;
        x_BB = 20; w_BB = 15;
        r_leftBB = Rect(-x_BB, 0, w_BB, w_BB,'base','center');
        r_rightBB = Rect(x_BB, 0, w_BB, w_BB,'base','center');
        layer_metalBB = dev.P_electrode.layer.metal1_BB;
        r_leftBB.layer = layer_metalBB;
        r_rightBB.layer = layer_metalBB;
        
        gtether_EC.addelement(r_elec1);
        gtether_EC.addelement(r_elec2);
        gtether_EC.addelement(r_elec_abs);
        gtether_EC.addelement(r_leftBB);
        gtether_EC.addelement(r_rightBB);
    end
    gtether_EC.rotate(theta);
    gtether_EC.translate(arc.p_end + L_halfTether*[-sind(theta);cosd(theta)]);
    devAll.addelement(gtether_EC);
end
% end coupler
v1 = [w/2;0]; v4 = [-w/2;0]; v2 = [w_end/2;L_endCoupler]; v3 = [-w_end/2;L_endCoupler];
p_EC = Polygon({v1,v2,v3,v4});
p_EC.layer = dev.layer_neg;
w_cure = 3;
l_cure_extra = 10; l_buff = 5;  % extra length for cure to include tolerance and cover the coupler holder
v1 = [w/2+w_cure;-l_cure_extra]; v4 = [-w/2-w_cure;-l_cure_extra];
v2 = [w_end/2+w_cure;L_endCoupler+l_buff]; v3 = [-w_end/2-w_cure;L_endCoupler+l_buff];
p_EC_cure = Polygon({v1,v2,v3,v4});
p_EC_cure.layer = dev.layer_cure;
% undercut mask
w_undercut = 5;
v1 = [w/2+w_undercut;-l_cure_extra]; v4 = [-w/2-w_undercut;-l_cure_extra];
v2 = [w_end/2+w_undercut;L_endCoupler+l_buff]; v3 = [-w_end/2-w_undercut;L_endCoupler+l_buff];
p_EC_undercut = Polygon({v1,v2,v3,v4});
p_EC_undercut.layer = dev.P_undercutMsk.layer;
g_EC = gpack.Group(0,0,{p_EC,p_EC_cure,p_EC_undercut});
g_EC.rotate(theta);
if abs(theta) > theta_genExtraTether
    r_disp = arc.p_end + 2*L_halfTether*[-sind(theta);cosd(theta)];
else
    r_disp = arc.p_end;
end
g_EC.translate(r_disp);
devAll.addelement(g_EC);

%% translate the whole device
if abs(theta) > theta_genExtraTether
    r_EC = arc.p_end + (2*L_halfTether+L_endCoupler)*[-sind(theta);cosd(theta)];
else
    r_EC = arc.p_end + (L_endCoupler)*[-sind(theta);cosd(theta)];
end
devAll.translate(-r_EC);
devAll.rotate(-theta);
dev.elements{dev.ind_expoBox}.translate([0;-7]);
dev.elements{dev.ind_expoBox}.rotate(theta);


%% big beam electrode
if isfield(param, 'isGenMetalBB') && param.isGenMetalBB
    layer_metalBB = dev.P_electrode.layer.metal1_BB;
    y_NBCenter = dev.y;
    l_NB = dev.l_NB_unscaled * dev.scale_NB;
    if theta >= 0
        sgn = 1;
    else
        sgn = -1;
    end    
    port1 = Points2D([[sgn*16.5;y_NBCenter + l_NB/2],...
        [sgn*(16.5 + 16.5); y_NBCenter + l_NB/2]]);
    port2 = Points2D([[sgn*20;y_NBCenter - l_NB/2],...
        [sgn*20; y_NBCenter - l_NB/2 - 16.5]]);    
    port1.translate(-r_EC);
    port1.rotate(-theta);   
    port2.translate(-r_EC);
    port2.rotate(-theta);
    
    y_pad = -300;
    w_pad = 250;   
    w_wire = 10;
    d_line2pad = 10 + w_wire/2;
    % wiring 2, from bottom electrode to bonding pad
    w2p1 = port2.points(:,1); w2p2 = port2.points(:,2);
    w2p3 = [0; y_pad];
    w2 = Wire({w2p1,w2p2,w2p3},w_wire);
    w2.layer = layer_metalBB;
    % wiring 1, no bonding pad
    w1p1 = port1.points(:,1); w1p2 = port1.points(:,2);
    y_line = y_pad + w_pad/2 + d_line2pad;
    x_line = w1p2(1) - w2p2(1); x2_line = sgn*(w_pad/2 + d_line2pad);
    y2_line = y_pad - w_pad/2 - d_line2pad;
    w1p3 = [x_line; y_line]; w1p4 = [x2_line; y_line]; w1p5 = [x2_line; y2_line];
    w1 = Wire({w1p1,w1p2,w1p3,w1p4,w1p5},w_wire);
    w1.layer = layer_metalBB;
    rect_bus1 = Rect(0, y2_line, 75*4*3, w_wire,'base','center');
    rect_bus1.layer = layer_metalBB;
    g_pad = genBondPad(0, y_pad, w_pad);
    g_BB = gpack.Group(0,0,{w1,w2,g_pad, rect_bus1});
    
    %
    devAll_BB = gpack.Group(0,0,{});
    devAll_BB.addelement(devAll);
    devAll_BB.addelement(g_BB);
    devAll = devAll_BB;
end

end



function ind = genLNNB_DS_EC_ind(ii)
global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        config.mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = -5-8+3;
dy = -12 - 6 - 3;
ind.rotate(-90);
ind.translate(dx,dy);

end