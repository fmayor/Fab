function tpnc = genTPNCDev(param)
% generate single TPNC device
% WTJ, 20190226

if param.flag_def
    tpnc = TPNC_2;
else
    tpnc = TPNC;
end
tpnc.P_PNC.w2 = 0.075;
tpnc.setN_PNC(3);
tpnc.P_PNC.N_mirror = 6;
tpnc.L_Bender = param.L_Bender;
tpnc.g_PS2PNCL = param.g_PS2PNCL;
tpnc.g_PS2PNCR = param.g_PS2PNCR;
tpnc.dy_PNC=8;
w_def = param.w_def;

d_w_def = w_def - 1.2;
ws_def = [1.2, 1.1, 1.3, 1.21, 1.25, 1.3]+d_w_def;
ws_def = num2cell(ws_def);
[tpnc.Ps_PNC_def.w1w_def] = ws_def{:};

% match layer names for Felix's convention
tpnc.layer_HD = 'M1_Neg';
tpnc.layer_MT = 'metal';
tpnc.layer_cure = 'cure';
tpnc.layer_RM = 'M2_undercutMsk';
tpnc.layer_MT_PL = 'metal_BB';

tpnc.run;
tpnc.addelement(genLNNB_DS_EC_ind(param.ind));


function g_ind = genLNNB_DS_EC_ind(ii)
% global config;
[ind_pos, ind_neg] = genStringPolygon(sprintf('%d',ii),...
        4);
ind_pos.layer = tpnc.layer_MT_PL;
ind_neg.layer = [tpnc.layer_MT_PL '_Neg'];
g_ind = gpack.Group(0,0,{ind_pos, ind_neg});
dx = 50;
dy = -250;
g_ind.rotate(-90);
g_ind.translate(dx,dy);

end



end