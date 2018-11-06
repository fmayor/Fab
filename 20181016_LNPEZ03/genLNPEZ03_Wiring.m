% genLNNB14 DC wiring


devMap = gpack.Group(0,0,{});

%% generate global marks
dx = 3500;
dy = 2000;
layers.pos = 'M1_Neg'; layers.neg = 'M1_Pos'; layers.M2 = 'M2_Neg';
GMs = genGlobalMarks(dx, dy, 3, 200, layers);
devMap.addelement(GMs);

%%
n_wiredDev = 13;
d_wiredDev = 62.5 * 4;

params = genParamStructs('d_nb2refl', [1.85],...
    'scale_NB', [0.99:0.06:1.11],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'isGenPS',false,'w2',0.05,...
    'isGenElectrode',false,'ElectrodeDy', 0.15,'d_metal2defect', 10);
dx = 62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
x_tot = dx * nx;
fprintf('Generating NB wiring...\n');
tic;
LNNBEC_wire = genPatternAry(params, xs, ys, @genSingleLNNBWire);
toc;


% add external bonding pad
% port 1, channel with onchip wirebond
x_bp1 = 3000/2 + 250 + 840;
x_bp2 = x_bp1 + 1500 - 840;
y_bp1 = -350;
bp1 = genBondPad(x_bp1, y_bp1);
bp2 = genBondPad(x_bp2, y_bp1, 250);
w = 10;
p_bp1ToBp2 = Wire({[x_bp1;y_bp1],[x_bp2;y_bp1]},w);
p_bp1ToBp2.layer = 'metal_BB';

% port 2 (channel without onchip wirebond)
x_port2 = 1600 + 840;
y_port2 = -550;
x_bondPad = x_port2 + 1650 - 840;
y_bondPad = y_port2 - 400;
w_buffer = 50;
r_buffer = Rect(x_port2, y_port2, w_buffer, w_buffer, 'base','center');
r_buffer.layer = 'metal_BB';
p_port2BondPad = Wire({[x_port2; y_port2],[x_bondPad; y_bondPad]},w);
p_port2BondPad.layer = 'metal_BB';
bp_port2 = genBondPad(x_bondPad, y_bondPad, 250);

g_metal_BB = gpack.Group(0,0,{r_buffer,bp1,bp2,bp_port2, p_bp1ToBp2, p_port2BondPad});
LNNBEC_wire.addelement(g_metal_BB);

LNNBEC_wire.translate(0,1700);


LNNBEC_wire2 = LNNBEC_wire.copy();
LNNBEC_wire2.rotate(180);

%%
devMap.addelement(LNNBEC_wire);
devMap.addelement(LNNBEC_wire2);
devMap.todxf('LNPEZ03_Wiring_1801018');

%% single device wiring
function dev = genSingleLNNBWire(dumb)
% for LNNB14
% WTJ, 20180719
ind = dumb.ind;
if mod(ind,4)==3
        
    [g_metal1, g_metal2] = genRightWire();
    [g_metal1_L, g_metal2_L] = genRightWire();
    g_metal1_L.mirror([0;0],[0;1]);
    g_metal2_L.mirror([0;0],[0;1]);
    
    dev = gpack.Group(0,0,{g_metal1, g_metal2, g_metal1_L,g_metal2_L});
else
    dev = gpack.Group(0,0,{});
end

end
function [g_metal1, g_metal2] = genRightWire()

    % metal type 1
    % top pad -> vertex 1 -> 
    w = 5;    % wire width
    x_pad = 25; y_topPad = -19;  y_botPad = -50;
    v0 = [x_pad; y_topPad]; v0Bot = [x_pad; y_botPad];
    v1 = [40; -19];  v2 = [40; -180];  v3 = [125; -180];  v4 = [125; -550];
    v5 = [0; -550];
    x_bondPad = 0; y_bondPad = -350;
    w_pad = 10;
    r_topPad = Rect(x_pad, y_topPad, w_pad, w_pad, 'base','center');
    r_botPad = Rect(x_pad, y_botPad, w_pad, w_pad, 'base','center');
    % dunno why this wire width is actually half as given
    p_topPad2port = Wire({v0,v1,v2,v3,v4, v5}, 2*w);
    p_botPad2bondPad = Wire({v0Bot, [x_bondPad; y_bondPad]},w);
    w_bondPad = 190;
    w_bondPadWindow = 180;
    r_bondPad = Rect(x_bondPad, y_bondPad, w_bondPad, w_bondPad, 'base','center');
    r_bondPadWin = Rect(x_bondPad, y_bondPad, w_bondPadWindow, w_bondPadWindow, 'base','center');
    r_bondPadWin.layer = 'M2_undercutMsk';
    
    g_metal1_metal = gpack.Group(0,0,{r_topPad, r_botPad, p_topPad2port,p_botPad2bondPad, r_bondPad});
    g_metal1_metal.layer = 'metal1_BB';
    g_metal1 = gpack.Group(0,0,{g_metal1_metal, r_bondPadWin});
    
    
    % metal type 2
    w = 5;    % wire width
    v0_p = [25.5; -42]; v0_n = [25.5; -56];
    v1_p = [40; -44.5];
    w_pad = 5;
    r_topPad = Rect(v0_p(1), v0_p(2), w_pad, w_pad, 'base','center');
    r_botPad = Rect(v0_n(1), v0_n(2), w_pad, w_pad, 'base','center');
    % dunno why this wire width is actually half as given
    p_topPad2port = Wire({v0_p,v1_p,v2,v3,v4, v5}, 2*w);
    p_botPad2bondPad = Wire({v0_n, [x_bondPad; y_bondPad]},w);
    r_bondPad = Rect(x_bondPad, y_bondPad, w_bondPad, w_bondPad, 'base','center');
    r_bondPadWin = Rect(x_bondPad, y_bondPad, w_bondPadWindow, w_bondPadWindow, 'base','center');
    r_bondPadWin.layer = 'M2_undercutMsk';
    
    g_metal2_metal = gpack.Group(0,0,{r_topPad, r_botPad, p_topPad2port,p_botPad2bondPad, r_bondPad});
    g_metal2_metal.layer = 'metal2_BB';
    g_metal2 = gpack.Group(0,0,{g_metal2_metal, r_bondPadWin});
end

function dev = genBondPad(x,y, w_bondPad)

if nargin < 3
    w_bondPad = 190;
end
    w_bondPadWindow = w_bondPad - 10;
    r_bondPad = Rect(0, 0, w_bondPad, w_bondPad, 'base','center');
    r_bondPad.layer = 'metal_BB';
    r_bondPadWin = Rect(0, 0, w_bondPadWindow, w_bondPadWindow, 'base','center');
    r_bondPadWin.layer = 'M2_undercutMsk';
    dev = gpack.Group(x,y,{r_bondPad,r_bondPadWin});
end

