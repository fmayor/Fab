% generate TPNC devices
% WTJ, 20190226


%%
% I have 5000 by 2200 um, could fit in ~ 10 by 3 TPNC devices
dx = 400;
dy = 600;

nx = 10;
ny = 3;
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = linspace(-dy*(ny-1)/2,dy*(ny-1)/2, ny);
%meshgrid is done within genPatternAry
params = genParamStructs('L_Bender',[5, 15, 25], 'g_PS2PNC', [0.2, 0.3],...
    'w_def', linspace(1.1, 1.2, 5));

tic;
TPNC_matrix = genPatternAry(params, xs, ys, @genTPNCDev);
toc;



%%
fname = 'TPNC02_20190307';
TPNC_matrix.todxf(fname);


