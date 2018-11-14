

function dev = genLNNB_DS_EC(param)
% generate single LNNB double-sided edge coupling device
global config;
w_end = param.w_end;
d_nb2refl = param.d_nb2refl;
scale_NB = param.scale_NB;
ind = param.ind;
tether_w2 = param.tether_w2;
tether_L = param.tether_L;
tether_w_tether = param.tether_w_tether;


% mirror cell from FBH GA 180707
P_mirror = struct('a', 557.6e-3, 'w', 1083e-3, 'amp', 81.9e-3, 'isFBH',true,...
    'hx', 330e-3, 'hy', 859e-3);
% best cost:
% P_defect = struct('a', 480.7e-3, 'w', 1344e-3, 'amp', 277.5e-3, 'isFBH',true,...
%     'hx', 352e-3, 'hy', 548.8e-3);
% ind = 207 from GA 180723_2
P_defect = struct('a', 450e-3, 'w', 1170e-3, 'amp', 39e-3, 'isFBH',true,...
    'hx', 334e-3, 'hy', 811e-3);

dev = nanobeam_DS_EC_SR;
dev.P_defect = P_defect;
dev.P_mirror = P_mirror;
dev.N_trans = 12;
dev.N_mirror = 17;
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
dev.layer_anchor = 'M1_LD';
dev.w_end = w_end;

dev.P_undercutMsk.isGen = 1;
dev.P_undercutMsk.layer = 'M2_undercutMsk';
dev.P_undercutMsk.d = 6.5;
% remove old buffer, new buffer will be added in addBendedEdgeCoupler
dev.P_undercutMsk.h = 0;    

dev.P_tether.w2 = tether_w2;
dev.P_tether.L = tether_L;
dev.P_tether.w_tether = tether_w_tether;
dev.P_tether.w_tetherArm = 17.5;    % tether length
% sharp
dev.P_sharp.isGen = true;
dev.P_sharp.layer = 'M1_Pos_holes';
dev.P_sharp.dutycyc = 0.9;
dev.P_sharp.period = 0.3;
dev.P_sharp.theta = 15;


% if mod(ind,4)==3
    dev.P_electrode.isGen = true;
    dev.P_electrode.dx = 22;
    dev.P_electrode.dy = 0.15;
    dev.P_electrode.d_metal2defect = 11;
    dev.P_electrode.layer.metal1 = 'metal';
    dev.P_electrode.layer.metal2 = 'trash';
    dev.P_electrode.layer.metal2_BB = 'trash';
    dev.P_electrode.layer.metal1_BB = 'metal_BB';
    dev.P_electrode.layer.metal1_abs = 'metal_abs';
% end
if param.isGenPS
    dev.P_1DPS.isGen = true;
    dev.P_1DPS.w1x = 0.6;
    dev.P_1DPS.w1y = 0.65;
    dev.P_1DPS.w2 = param.w2;
end

dev.scale_NB = scale_NB;
dev.maker = @FBHMaker_F;
dev.P_coupler.maker = @FBHMaker_F;

% apply correction
correct_LNNB16(dev);
% do not run, wait for adding bended edge coupler
% dev.run(config.isfast);
% % add index
% g_ind = genLNNB_DS_EC_ind(ind);
% dev.addelement(g_ind);


% add bended coupler
dev = addBendedEdgeCoupler(dev,param);

end

