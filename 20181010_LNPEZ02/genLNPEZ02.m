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
%% generate LNNB22
% 74 devices in total
params = genParamStructs('d_nb2refl', [1.75, 1.8],...
    'scale_NB', [1.04:0.002:1.12],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'isGenPS',false,'w2',0.05);
% phonon shield parameters
for ii = 1:length(params)
    if mod(ii,3) ~= 0
        params(ii).isGenPS = true;
        if mod(ii,3) == 2
            params(ii).w2 = 0.065;   
%             params(ii).w2 = 0.04;
        end
    end
end
% thetas = [-flip(linspace(0,90,7))-15/2, linspace(0,90,7)];
% thetas = thetas(2:end);
% thetas = [thetas -thetas thetas -thetas thetas -thetas];
thetas = [0,30,40,45,50,60,90,75,55,47.5,45, 42.5, 35,15];
thetas = [thetas -thetas thetas -thetas thetas -thetas];
thetas = num2cell(thetas(1:length(params)));
[params.theta] = thetas{:};
dx = 75;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
x_tot = dx * nx;
fprintf('Generating NBs...\n');
tic;
LNNBEC_Row = genPatternAry(params, xs, ys, @genLNNB_DS_EC);
toc;
LNNBEC_Row.translate(0, 1700);
devMap.addelement(LNNBEC_Row);
% LNNBEC_Row_bot = LNNBEC_Row.copy(); LNNBEC_Row_bot.rotate(180);
% devMap.addelement(LNNBEC_Row_bot);

%% generate PEZZS devices
% sweeping zigzag and mid coupler parameters
params1 = genParamStructs('w_mid', [0.3,0.4,0.5,0.7],...
    'w_end', 0.26, 'L', [5:5:20],...
    'L2',[5:5:10],'g', [0.3], 'g_mid', [0, 0.3]);
% % control devices 1
% params1_c = genParamStructs('w_mid', [0.3 : 0.1 : 0.8],...
%     'w_end', 0.26, 'L', [5],...
%     'L2',[1],'g', [0.3], 'g_mid', [0, 0.3]);
% % sweeping g_mid
params2 = genParamStructs('w_mid', 0.4,...
    'w_end', 0.26, 'L', 5,...
    'L2',2,'g', [0.3], 'g_mid', [0:0.05:0.5]);
% % control devices 2
% params2_c = genParamStructs('w_mid', [0.3 : 0.1 : 0.8],...
%     'w_end', 0.26, 'L', [5],...
%     'L2',[1],'g', [0.3], 'g_mid', [0, 0.3]);
params = [params2 params1 params2];
dx = 62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
fprintf('Generating PEZZSs...\n');
tic;
PEZZS_Row = genPatternAry(params, xs, ys, @genPEZZSsym);
toc;
PEZZS_Row.translate(0, 1700);
PEZZS_Row.rotate(180);
devMap.addelement(PEZZS_Row);

%% generate Dev for SEM
param = struct('scale_NB', 1,...
    'd_nb2refl', 1.8,'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[4],'tether_w_tether', [1],'isGenPS',true,'w2',0.05);
param.ind = 0;
dev4SEM = genLNNB_DS_EC(param);
g_deb4SEM = gpack.Group(0,0,{dev4SEM});
g_deb4SEM_bot = g_deb4SEM.copy();
g_deb4SEM.translate(0,1000);
g_deb4SEM_bot.translate(0,-1000);
devMap.addelement(g_deb4SEM);
devMap.addelement(g_deb4SEM_bot);



%% generate zigzags

% params = genParamStructs('w', [0.300, 0.400], 'g', [0.3, 0.4],...
%                     'L2',[10, 20, 30], 'w_joint', [0.03, 0.050]);
% 
% %
% nx = length(params);
% dx = 62.5;
% xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
% ys = [0];
% fprintf('Generating zigzags...\n');
% zigzags = genPatternAry(params, xs, ys, @genDoubleZZBender_180628);
% zigzags2 = zigzags.copy();
% zigzags.translate(0,100);
% zigzags2.translate(0,-100);
% dev = gpack.Group(0,0,{zigzags,zigzags2});
% devMap.addelement(dev);

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

fname = 'LNPEZ_20181010';
if config.isfast
    fname = [fname '_fast'];
end

if isnoGC
    fname = [fname '_noGC'];
end
devMap.todxf(fname);

%% single device functions
