function dev = genLNZNB_DS_EC_ZZ_X(param)
% generate single LN zipper NB attached to zigzag double-sided edge coupling device 
% 
% WTJ, FMM, 20181128


global config;
w_end = param.w_end;
d_nb2refl = param.d_nb2refl;
scale_NB = param.scale_NB;
ind = param.ind;
tether_w2 = param.tether_w2;
tether_L = param.tether_L;
tether_w_tether = param.tether_w_tether;
gap = param.gap;


theta = 22/180*pi;
theta_out = 12/180*pi;

% mirror cell from model_FBHZUC_181103.mph
P_mirror = struct('a', 650e-3, 'w', 1351e-3, 'amp', 97.1e-3, 'isFBH',true,...
    'hx', 480.8e-3, 'hy', 1102.6e-3,'theta',theta, 'theta_out',theta_out,...
    'gap', 150e-3);
% ind = 12 from GA 181103
% P_defect = struct('a', 543.9e-3, 'w', 1460e-3, 'amp', 42.7e-3, 'isFBH',true,...
%     'hx', 289.3e-3, 'hy', 934.6e-3,'theta',theta, 'theta_out',theta_out,...
%     'gap', 150e-3);
% ind = 68 from GA 181103
P_defect = struct('a', 546.9e-3, 'w', 1512e-3, 'amp', 16.6e-3, 'isFBH',true,...
    'hx', 381e-3, 'hy', 843e-3,'theta',theta, 'theta_out',theta_out,...
    'gap', 150e-3);
P_mirror.gap = gap;
P_defect.gap = gap;
if param.theta == 0
    dev = nanobeam_DS_EC_ZZ_X_probes_rotated;
else
    dev = nanobeam_DS_EC_ZZ_X_probes;
end
%dev = nanobeam_DS_EC_ZZ_X_probes;
dev.w=param.w;%zigzag parameters
dev.L=param.L;
dev.L2=param.L2;
dev.g=param.g;
dev.g_metal=param.g_metal;
dev.P_defect = P_defect;
dev.P_mirror = P_mirror;
dev.N_trans = 12;
dev.N_mirror = 13+round(dev.L/P_mirror.a);
dev.d_tether2edge = 10+7;
dev.P_coupler.l_extra = 0;    % coupler end at buffer box edge
dev.P_tether.L = 3.5;
dev.P_tether.w_tether = 1;
dev.d_nb2edge = 7;
dev.d_nb2refl = d_nb2refl;
dev.w_buffer = 62.5;
dev.layer_pos = 'M1_Pos';
dev.layer_neg = 'M1_Neg';
dev.layer_wf = 'M1_Field';
dev.layer_holes = 'M1_Pos_holes';
dev.layer_medium = 'M1_MD';
dev.layer_metal_BB = 'metal_BB';
dev.w_end = w_end;
dev.P_undercutMsk.isGen = 1;
dev.P_undercutMsk.layer = 'M2_undercutMsk';
dev.P_undercutMsk.d = 4.5;
dev.h_link = 0.4;
dev.w_link = 0.05;
dev.P_tether.w2 = tether_w2;
dev.P_tether.L = tether_L;
dev.P_tether.w_tether = tether_w_tether;
dev.P_tether.w_tetherArm = 2*(d_nb2refl+P_mirror.gap/2+P_mirror.w/2+P_mirror.amp/2+dev.w_link);%10;
if mod(param.ind,2)==0
    dev.isSB = true;
    dev.isBB = true;
else
    dev.isSB = false;
    dev.isBB = false;
end

dev.layer_anchor = 'M1_LD';
if param.isGenPS
    dev.P_1DPS.isGen = true;
    dev.P_1DPS.w1x = 0.6;
    dev.P_1DPS.w1y = 0.65;
    dev.P_1DPS.w2 = param.w2;
end

dev.scale_NB = scale_NB;
dev.maker = @FBHZMaker_F;
dev.P_coupler.maker = @FBHZMaker_F;

% apply correction
correct_LNZNB(dev);
dev.run(config.isfast);


%add bended coupler to align ZZ arm with crystal X
dev=addBendedEdgeCoupler(dev,param);

% make M1_LD within M1_Field
dev.elements{1}.elements{dev.elements{1}.ind_expoBox}.translate([13.5-7;0]);

%make small electrodes
% port_RL=dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{2}.port1;
% port_RL = port_RL+[dev.x+dev.elements{1}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{2}.x;dev.y+dev.elements{1}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{2}.y];
% port_RR=dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{3}.port1;
% port_RR = port_RR+[dev.x+dev.elements{1}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.x+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{3}.x;dev.y+dev.elements{1}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.y+dev.elements{1}.elements{dev.elements{1}.indZZ1}.elements{2}.elements{3}.y];
% 
% %right zigzag electrodes
% param_elec.dy=6;
% param_elec.w_pad=4;
% param_elec.zz_w=param.w;
% param_elec.port_L=port_RL;
% param_elec.port_R=port_RR;
% [small_R_elec,port_B_RL,port_B_RR]=genSmallAnchorElectrodes(param_elec);
% %small_R_elec.translate([-2*dx;0]);
% dev.addelement(small_R_elec);




% add index
if param.theta == 0
    g_ind = genLNNB_DS_EC_ind_rotated(ind);

else
    g_ind = genLNNB_DS_EC_ind(ind);

end
%g_ind = genLNNB_DS_EC_ind(ind);
dev.addelement(g_ind);

%add UC mask
if param.theta == 0
    UC_mask=Rect(0,10,260,20,'base','center');

else
    UC_mask=Rect(55.5,-34,20,260,'base','center');

end
UC_mask.layer='M2_undercutMsk';
dev.addelement(UC_mask);
end

function ind = genLNNB_DS_EC_ind(ii)
global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        config.mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = 5+8+5+5;
dy = -12 - 20-8+15+2+7.5;
ind.rotate(180);
ind.translate(dx,dy);

end

function ind = genLNNB_DS_EC_ind_rotated(ii)
global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        config.mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = 5+8+5+5+10;
dy = -12 - 20-8+15+7;
ind.rotate(-90);
ind.translate(dx,dy);

end