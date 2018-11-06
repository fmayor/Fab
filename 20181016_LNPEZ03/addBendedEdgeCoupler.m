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
    ind = param.ind;
    g_ind = genLNNB_DS_EC_ind(ind);
    dev.addelement(g_ind);
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
ind = param.ind;
g_ind = genLNNB_DS_EC_ind(ind);
dev.addelement(g_ind);
devAll.addelement(dev);
%% generate arc
arc = Arc_zero([0;0], w, R, theta);
arc.layer = dev.layer_neg;
arc_cure = Arc_zero([0;0], w + 6, R, theta);
arc_cure.layer = dev.layer_cure;
devAll.addelement(arc);
devAll.addelement(arc_cure);

%% generate tether and end-coupler
% tether
theta_genExtraTether = 30;
if abs(theta) > theta_genExtraTether
gtether = dev.elements{dev.ind_gtether};
gtether_EC = gtether.copy();
gtether_EC.y = 0;
gtether_EC.rotate(theta);
gtether_EC.translate(arc.p_end + L_halfTether*[-sind(theta);cosd(theta)]);
devAll.addelement(gtether_EC);
end
% end coupler
v1 = [w/2;0]; v4 = [-w/2;0]; v2 = [w_end/2;L_endCoupler]; v3 = [-w_end/2;L_endCoupler];
p_EC = Polygon({v1,v2,v3,v4});
p_EC.layer = dev.layer_neg;
w_cure = 3;
v1 = [w/2+w_cure;0]; v4 = [-w/2-w_cure;0];
v2 = [w_end/2+w_cure;L_endCoupler]; v3 = [-w_end/2-w_cure;L_endCoupler];
p_EC_cure = Polygon({v1,v2,v3,v4});
p_EC_cure.layer = dev.layer_cure;
g_EC = gpack.Group(0,0,{p_EC,p_EC_cure});
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

end



function ind = genLNNB_DS_EC_ind(ii)
global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        config.mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = -5-8+3;
dy = -12 - 6;
ind.rotate(-90);
ind.translate(dx,dy);

end