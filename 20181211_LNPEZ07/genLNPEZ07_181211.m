% layout LNPEZ07



devMap = gpack.Group(0,0,{});

%% generate global marks
dx = 3000;
dy = 2000; %2000
layers.pos = 'M1_Neg'; layers.neg = 'M1_Pos'; layers.M2 = 'M2_Neg';
GMs = genGlobalMarks(dx, dy, 3, 200, layers);
devMap.addelement(GMs);

CM = genGlobalMark(0,0,'GC');
devMap.addelement(CM);
%% setup config
global config;
% config.isfast = true;
config.isfast = false;
config.mkrPixelSize = 0.5;
% config.LNNB.dx = 62.5;  % separation by one writefield
isnoGC = true;
%% generate top row
 [LNPEZ_ZNB_Z_Row, L_NBZZ] = genLNPEZ_ZNB_X_toprow_181211();
 LNPEZ_ZNB_Z_Row.translate(0,1700);
 devMap.addelement(LNPEZ_ZNB_Z_Row);

%% generate bot row
 [LNPEZ_ZNB_Z_Row, L_NBZZ2] = genLNPEZ_ZNB_X_botrow_181211();
 LNPEZ_ZNB_Z_Row.translate(0,1700);
 LNPEZ_ZNB_Z_Row.rotate(180);
 devMap.addelement(LNPEZ_ZNB_Z_Row);


%% manually add BB wire
w_wire = 10;
w_pad = 300;
% 
% %%bot row
% 
% %bonding pad from chip to PCB
%right
v_pad = [3000;-840];
w1 = Wire({[2770;-636],v_pad},w_wire);
v_pad2 = [3000; -460];
w2 = Wire({[2750;-420],v_pad2},w_wire);
g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
%left
v_pad3 = [-3100;-926];
w3 = Wire({[-2930;-780],v_pad3},w_wire);
v_pad4 = [-3070; -575];
w4 = Wire({[-2858;-641],v_pad4},w_wire);
g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);
w5= Wire({[-2800;-390],[-3271;-390],[-3271;-937],[-3028;-937]},16);


g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4,w5});
g_wire_BB.layer = 'metal_BB';
devMap.addelement(g_pad1);
devMap.addelement(g_pad2);
devMap.addelement(g_pad3);
devMap.addelement(g_pad4);
devMap.addelement(g_wire_BB);
%% top row

%right
v_pad = [3100;926];
w1 = Wire({[2979;781],v_pad},w_wire);
v_pad2 = [3070; 575];
w2 = Wire({[2847;642],v_pad2},w_wire);
g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
w5= Wire({[2800;390],[3271;390],[3271;937],[3028;937]},16);
%left
v_pad3 = [-3000;460];
w3 = Wire({[-2754;421],v_pad3},w_wire);
v_pad4 = [-3000; 840];
w4 = Wire({[-2761; 637],v_pad4},w_wire);
g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);
 
g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4,w5});
g_wire_BB.layer = 'metal_BB';
devMap.addelement(g_pad1);
devMap.addelement(g_pad2);
devMap.addelement(g_pad3);
devMap.addelement(g_pad4);
devMap.addelement(g_wire_BB);

%% add dicing mark
w_DM = 50;
dx = 3500;
dy = 1700 + 7 + w_DM/2;
dm1 = Rect(0,0,w_DM, w_DM,'base','center');
% dm1 = genTurkey_1801106(40);
dm1.layer = 'M1_LD';   % LD: low dose
dm2 = dm1.copy(); dm3 = dm1.copy(); dm4 = dm1.copy();
dm1.translate(-dx, dy);
dm2.translate(dx, dy);
dm3.translate(dx, -dy);
dm4.translate(-dx, -dy);
devMap.addelement(dm1);
devMap.addelement(dm2);
devMap.addelement(dm3);
devMap.addelement(dm4);

% add dektak marks
dx = 2140;
dy = 1044;
% dm1 = Rect(0,0,w_DM, w_DM,'base','center');
tr1 = genXmastree_1801203(35);
tr2 = genXmastree_1801203(35);
tr3 = genXmastree_1801203(35);
tr2.translate([11.5;0]);
tr3.translate([-11.5;0]);
dm1=gpack.Group(0,0,{tr1,tr2,tr3});
dm1.layer = 'M1_LD';
dm2 = dm1.copy(); dm3 = dm1.copy(); dm4 = dm1.copy();
dm1.translate(-dx, dy);
dm2.translate(dx, dy);
dm3.translate(dx, -dy+16);
dm4.translate(-dx, -dy+16);
devMap.addelement(dm1);
devMap.addelement(dm2);
devMap.addelement(dm3);
devMap.addelement(dm4);

%% to dxf

fname = 'LNPEZ07_181211_M1';
if config.isfast
    fname = [fname '_fast'];
end

if isnoGC
    fname = [fname '_noGC'];
end
devMap.todxf(fname);

%% single device functions


function GM = genGlobalMark(dx,dy,name)
layers = struct('pos','M1_Neg','neg','M1_Pos','M2','M2_Neg');
GM = NestedAlignmentMark([0, 0], layers.pos,...
    0.4,2,4,60);
mkrPixelSize = 0.7;
[GMInd_pos, GMInd_neg] = genStringPolygon(name, mkrPixelSize);
GMInd_pos.layer = layers.pos; GMInd_neg.layer = layers.neg;
GMInd_pos.rotate(-90); GMInd_neg.rotate(-90);
dx_ind = -35; dy_ind = -5;
GMInd_pos.translate(dx_ind,dy_ind); GMInd_neg.translate(dx_ind,dy_ind);
GM.addelement(GMInd_pos); GM.addelement(GMInd_neg);
GM.translate(dx,dy);
% photolitho XeF2 mask
w_msk = 60*2+60;
GM_msk = Rect(0,0,w_msk,w_msk,'base','center');
GM_msk.layer = layers.M2;
GM.addelement(GM_msk);
end
