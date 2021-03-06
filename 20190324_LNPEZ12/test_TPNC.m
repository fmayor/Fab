% test TPNC
% WTJ, 20190226

%% test PNC



%% test Bender2PS
bdr2ps = Bender2PS_90deg;
bdr2ps.w_tip = 0.05;
bdr2ps.l_tip = 1.2;
% bdr2ps.cfg_wiring = 'B';
bdr2ps.run;
bdr2ps.todxf('tmp_bdr2ps');

%%

tpnc = TPNC;
tpnc.P_PNC.w2 = 0.075;
tpnc.L_Bender = 15;
tpnc.setN_PNC(4);
ws_def = [1.2, 1.1, 1.2, 1.3];%[1.2, 1.1, 1.15, 1.2, 1.25, 1.3];
ws_def = num2cell(ws_def);
[tpnc.Ps_PNC_def.w1w_def] = ws_def{:};

tpnc.run;
% add available region
tpnc.addelement(Rect(0,0, 5000, 2200,'base','center'));
tpnc.todxf('test_tpnc');

%% test single dev fun
param.L_Bender = 5;
param.ind = 30;
param.w_def = 1.22;
param.g_PS2PNC = 0.2;
tpnc = genTPNCDev(param);
tpnc.todxf('test_tpnc');
