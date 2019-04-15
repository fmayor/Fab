% layout LNPEZ12



devMap = gpack.Group(0,0,{});

%% generate global marks
dx = 3000;
dy = 2000; %2000
layers.pos = 'M1_Neg'; layers.neg = 'M1_Pos'; layers.M2 = 'M2_Neg';
GMs = genGlobalMarks(dx, dy, 3, 200, layers);
devMap.addelement(GMs);

CM = genGlobalMark(0,0,'GC'); %central mark
devMap.addelement(CM);
dy2=550;
GJ = genGlobalMark(-dx, dy2,'GJ');
GK = genGlobalMark(dx, dy2,'GK');
GL = genGlobalMark(dx, -dy2,'GL');
GM = genGlobalMark(-dx, -dy2,'GM');
devMap.addelement(Group(0,0,{GJ,GK,GL,GM}));
%% setup config
global config;
% config.isfast = true;
config.isfast = false;
config.mkrPixelSize = 0.5;
% config.LNNB.dx = 62.5;  % separation by one writefield
isnoGC = true;
offset=500; %second row offset from first row
%% generate top row
 [LNPEZ_ZNB_Z_Row, L_NBZZ] = genLNPEZ_ZNB_X_toprow_190324();
 LNPEZ_ZNB_Z_Row.translate(0,1700);
 devMap.addelement(LNPEZ_ZNB_Z_Row);
LNPEZ_ZNB_Z_Row2=LNPEZ_ZNB_Z_Row.copy();
 LNPEZ_ZNB_Z_Row2.translate(0,-offset);
 devMap.addelement(LNPEZ_ZNB_Z_Row2);
%% generate bot row
 [LNPEZ_ZNB_Z_Row, L_NBZZ2] = genLNPEZ_ZNB_X_botrow_190324();
 LNPEZ_ZNB_Z_Row.translate(0,1700);
 LNPEZ_ZNB_Z_Row.rotate(180);
 devMap.addelement(LNPEZ_ZNB_Z_Row);
 LNPEZ_ZNB_Z_Row2=LNPEZ_ZNB_Z_Row.copy();
 LNPEZ_ZNB_Z_Row2.translate(0,offset);
 devMap.addelement(LNPEZ_ZNB_Z_Row2);

 %% generate mid row
%  [CL_Row, L_ZZCL] = genLNPEZ_ZNB_X_midrow_190224(); %cantilever ZZ
%  CL_Row.translate(0,150);
%  devMap.addelement(CL_Row);
%% manually add BB wire
% w_wire = 10;
% w_pad = 300;
% % 
% % %%bot row
% % 
% % %bonding pad from chip to PCB
% %right
% v_pad = [3000;-840];
% w1 = Wire({[2720;-636],v_pad},w_wire);
% v_pad2 = [3000; -460];
% w2 = Wire({[2700;-420],v_pad2},w_wire);
% g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
% g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
% %left
% v_pad3 = [-3100;-926];
% w3 = Wire({[-2930;-780],v_pad3},w_wire);
% v_pad4 = [-3070; -575];
% w4 = Wire({[-2825;-641],v_pad4},w_wire);
% g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
% g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);
% w5= Wire({[-2800;-390],[-3271;-390],[-3271;-937],[-3028;-937]},16);
% 
% 
% g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4,w5});
% g_wire_BB.layer = 'metal_BB';
% devMap.addelement(g_pad1);
% devMap.addelement(g_pad2);
% devMap.addelement(g_pad3);
% devMap.addelement(g_pad4);
% devMap.addelement(g_wire_BB);
%% top row
% 
% %right
% v_pad = [3100;926];
% w1 = Wire({[2979;781],v_pad},w_wire);
% v_pad2 = [3070; 575];
% w2 = Wire({[2847;642],v_pad2},w_wire);
% g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
% g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
% w5= Wire({[2800;390],[3271;390],[3271;937],[3028;937]},16);
% %left
% v_pad3 = [-3000;460];
% w3 = Wire({[-2754;421],v_pad3},w_wire);
% v_pad4 = [-3000; 840];
% w4 = Wire({[-2761; 637],v_pad4},w_wire);
% g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
% g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);
%  
% g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4,w5});
% g_wire_BB.layer = 'metal_BB';
% devMap.addelement(g_pad1);
% devMap.addelement(g_pad2);
% devMap.addelement(g_pad3);
% devMap.addelement(g_pad4);
% devMap.addelement(g_wire_BB);

%% add dicing mark
w_DM = 50;
dx = 3500;
dy = 1700 + 7 + w_DM/2;
dm1 = Rect(0,0,w_DM, w_DM,'base','center');
% dm1 = genTurkey_1801106(40);
dm1.layer = 'M1_LD';   % LD: low dose
dm2 = dm1.copy(); dm3 = dm1.copy(); dm4 = dm1.copy(); 
dm5 = dm1.copy(); dm6 = dm1.copy(); dm7 = dm1.copy(); dm8 = dm1.copy();
dm1.translate(-dx, dy);
dm2.translate(dx, dy);
dm3.translate(dx, -dy);
dm4.translate(-dx, -dy);
dm5.translate(-dx, dy-offset);
dm6.translate(dx, dy-offset);
dm7.translate(dx, -dy+offset);
dm8.translate(-dx, -dy+offset);
devMap.addelement(dm1);
devMap.addelement(dm2);
devMap.addelement(dm3);
devMap.addelement(dm4);
devMap.addelement(dm5);
devMap.addelement(dm6);
devMap.addelement(dm7);
devMap.addelement(dm8);
% add dektak marks
dx = 2140;
dy = 1044;
% dm1 = Rect(0,0,w_DM, w_DM,'base','center');
bunny = genBunny_190225(35);
dm1=gpack.Group(0,0,{bunny});
dm1.layer = 'M1_LD';
dm1.rotate(90);
dm2 = dm1.copy(); dm3 = dm1.copy(); dm4 = dm1.copy();
dm1.translate(-dx, dy-21);
dm2.translate(dx, dy-20);
dm3.translate(dx, -dy+16-20);
dm4.translate(-dx, -dy+16-21);
devMap.addelement(dm1);
devMap.addelement(dm2);
devMap.addelement(dm3);
devMap.addelement(dm4);
bunny = genBunny_190225(1); %for fun
dm5=gpack.Group(0,0,{bunny});
dm5.layer = 'M1_Neg';
dm5.rotate(90);
dm5.translate(-dx+50, dy-21);
devMap.addelement(dm5);
%% to dxf

fname = 'LNPEZ12_190328_M3';
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
