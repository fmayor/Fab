
function dev = genLNZNB_DS_EC(param)
% generate single LN zipper NB double-sided edge coupling device
% 
% WTJ, FMM, 20181106


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

dev = nanobeam_DS_EC_SR;
dev.P_defect = P_defect;
dev.P_mirror = P_mirror;
dev.N_trans = 12;
dev.N_mirror = 15;
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
dev.w_end = w_end;
dev.P_undercutMsk.isGen = 1;
dev.P_undercutMsk.layer = 'M2_undercutMsk';
dev.P_undercutMsk.d = 4.5;
dev.P_tether.w2 = tether_w2;
dev.P_tether.L = tether_L;
dev.P_tether.w_tether = tether_w_tether;
dev.P_tether.w_tetherArm = 10;

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
% add index
g_ind = genLNNB_DS_EC_ind(ind);
dev.addelement(g_ind);

% make M1_LD within M1_Field
dev.elements{dev.ind_expoBox}.translate(0,-2);

end



function ind = genLNNB_DS_EC_ind(ii)
global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        config.mkrPixelSize);
ind_pos.layer = 'M1_Neg';
ind_neg.layer = 'M1_Pos_ind';
ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = -5-8;
dy = -12 - 20-8;
ind.rotate(-90);
ind.translate(dx,dy);

end