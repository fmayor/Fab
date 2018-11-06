% layout LNNB09 



devMap = gpack.Group(0,0,{});

%% generate global marks
dx = 3500;
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
%% generate LNNBZZ
% 42 devices in total
params = genParamStructs('d_nb2refl', [1.75, 1.85],...
    'scale_NB', [1.02:0.005:1.12],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false);
params_ZZ = genParamStructs('g_metal',[0.12],'g_zigzag',0.3,'w_zigzag',[0.3, 0.45],...
    'zigzag_tether', [0.05, 0.1], 'L2', [10,15,20,15]);
params_ZZ = [params_ZZ params_ZZ params_ZZ];
params = copyfield(params, params_ZZ);
% reorder params randomly, generated once and then saved, so that it's
% repeatable
% rnds = rand(1,length(params));
% [~,inds] = sort(rnds);
load('randind_181017.mat');
params = params(inds);
% gen BB once every four devices shield parameters
for ii = 1:length(params)
    if mod(ii,4) == 0
        params(ii).isGenBB = true;
    end
end
dx = 62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
L_NBZZ = dx*nx;
fprintf('Generating LNNBZZs...\n');
tic;
LNNBZZ_Row = genPatternAry(params, xs, ys, @test_LNNB_Zigzag);
toc;
% devMap.addelement(LNNBZZ_Row);
% LNNBEC_Row_bot = LNNBEC_Row.copy(); LNNBEC_Row_bot.rotate(180);
% devMap.addelement(LNNBEC_Row_bot);

%% generate PEZZS devices
% sweeping zigzag and mid coupler parameters
params1 = genParamStructs('w_mid', [0.3,0.4,0.5,0.6],...
    'w_end', 0.26, 'L', [10,15,20],...
    'L2',[10],'g', [0.3], 'g_mid', [0.15],'zz_w',[0.3,0.45],'g_metal',[0.1,0.15]);

params = [params1];
dx = 80; %62.5
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
fprintf('Generating PEZZSs...\n');
tic;
PEZZS_Row = genPatternAry(params, xs, ys, @genPEZZSsym_wired);
toc;
L_PEZ=dx*nx;
L_tot = L_PEZ + L_NBZZ;
PEZZS_Row.translate(-(L_tot/2 - L_PEZ/2), 1700);
LNNBZZ_Row.translate((L_tot/2 - L_NBZZ/2), 1700);
%PEZZS_Row.rotate(180);
topRow = gpack.Group(0,0,{PEZZS_Row, LNNBZZ_Row});
devMap.addelement(topRow);


%% manually add BB wire
w_wire = 7;
w1 = Wire({[530;762],[680;747]},w_wire);
w2 = Wire({[530;519.6],[681;519.6]},w_wire);

% bonding pad from chip to PCB
%right
v_pad = [3520;965];
w3 = Wire({[3223;747],v_pad},w_wire);
v_pad2 = [3520; 519.6];
w4 = Wire({[3223;519.6],v_pad2},w_wire);
w_pad = 300;
g_pad1 = genBondPad(v_pad(1),v_pad(2), w_pad);
g_pad2 = genBondPad(v_pad2(1),v_pad2(2), w_pad);
%left
v_pad3 = [-3616;762+223];
w5 = Wire({[-3316;762],v_pad3},w_wire);
v_pad4 = [-3616; 519.6];
w6 = Wire({[-3316;519.6],v_pad4},w_wire);
w_pad = 300;
g_pad3 = genBondPad(v_pad3(1),v_pad3(2), w_pad);
g_pad4 = genBondPad(v_pad4(1),v_pad4(2), w_pad);


g_wire_BB = gpack.Group(0,0,{w1,w2,w3,w4,w5,w6});
g_wire_BB.layer = 'metal_BB';
topRow.addelement(g_wire_BB);
topRow.addelement(g_pad1);
topRow.addelement(g_pad2);
topRow.addelement(g_pad3);
topRow.addelement(g_pad4);

%% bot row
botRow = topRow.copy();
botRow.rotate(180);
devMap.addelement(botRow);

%% add dicing mark
w_DM = 50;
dx = 3500;
dy = 1700 + 7 + w_DM/2;
dm1 = Rect(0,0,w_DM, w_DM,'base','center');
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
dx = 2000;
dy = 1000;
dm1 = Rect(0,0,w_DM, w_DM,'base','center');
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

fname = 'LNPEZ_20181018_M3';
if config.isfast
    fname = [fname '_fast'];
end

if isnoGC
    fname = [fname '_noGC'];
end
devMap.todxf(fname);

%% single device functions



