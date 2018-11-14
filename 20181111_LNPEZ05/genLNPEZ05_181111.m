% layout LNPEZ04



devMap = gpack.Group(0,0,{});

%% generate global marks
dx = 3000;
dy = 2000;
layers.pos = 'M1_Neg'; layers.neg = 'M1_Pos'; layers.M2 = 'M2_Neg';
GMs = genGlobalMarks(dx, dy, 3, 200, layers);
devMap.addelement(GMs);

%% setup config
global config;
% config.isfast = true;
config.isfast = false;
config.mkrPixelSize = 0.5;
% config.LNNB.dx = 62.5;  % separation by one writefield
isnoGC = true;
%% generate LNNB
%  devices in total
[LNZNB_ZZ_X_Row, L_NBZZ] = genLNZNB_row_181111();
LNZNB_ZZ_X_Row.translate(0,1700);
devMap.addelement(LNZNB_ZZ_X_Row);

%% generate PEZZS devices
%  [LNPEZ_ZNB_Z_Row, L_NBZZ2] = genLNPEZ_ZNB_Z_row_181111();
%  LNPEZ_ZNB_Z_Row.translate(0,1700);
%  LNPEZ_ZNB_Z_Row.rotate(180);
%  devMap.addelement(LNPEZ_ZNB_Z_Row);


%% manually add BB wire
% w_wire = 7;
% w_pad = 300;
% 
% %bonding pad from chip to PCB
% %right
% v_pad = [3000;-1100];
% w1 = Wire({[2627;-587],v_pad},w_wire);
% v_pad2 = [3000; -636];
% w2 = Wire({[2500;-330],v_pad2},w_wire);
% g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
% g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
% %left
% v_pad3 = [-3000;-1100];
% w3 = Wire({[-2627;-587],v_pad3},w_wire);
% v_pad4 = [-3000; -636];
% w4 = Wire({[-2500;-330],v_pad4},w_wire);
% g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
% g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);
%  
% g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4});
% g_wire_BB.layer = 'metal_BB';
% devMap.addelement(g_pad1);
% devMap.addelement(g_pad2);
% devMap.addelement(g_pad3);
% devMap.addelement(g_pad4);
% devMap.addelement(g_wire_BB);
%% bot row, 


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
dy = 1035;
% dm1 = Rect(0,0,w_DM, w_DM,'base','center');
dm1 = genTurkey_1801106(50);
dm1.layer = 'M1_LD';
dm2 = dm1.copy(); dm3 = dm1.copy(); dm4 = dm1.copy();
dm1.translate(-dx, dy);
dm2.translate(dx, dy);
dm3.translate(dx, -dy);
dm4.translate(-dx, -dy);
devMap.addelement(dm1);
devMap.addelement(dm2);
devMap.addelement(dm3);
devMap.addelement(dm4);

%% to dxf

fname = 'LNPEZ05_181113_M2';
if config.isfast
    fname = [fname '_fast'];
end

if isnoGC
    fname = [fname '_noGC'];
end
devMap.todxf(fname);

%% single device functions



