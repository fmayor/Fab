% generate TPNC devices
% WTJ,FMM, 20190324


%%
% I have 5000 by 2200 um, could fit in ~ 12 by 2 TPNC devices
dx = 400;
dy = 600;

nx = 12;
ny = 2;
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = linspace(-dy*(ny-1)/2,dy*(ny-1)/2, ny);
%meshgrid is done within genPatternAry
params = genParamStructs('w_def', [1.12,1.16],'L_Bender',[15,20,25],...
'g_PS2PNCL',[0,0,0,0]); 
%make different gaps for left and right tips
g_PS2PNCL=[0,0.25,0.35,0.45,0.05,0.25,0.35,0.45,0.1,0.25,0.35,0.45];
g_PS2PNCL=[g_PS2PNCL,repmat([0.25,0.3,0.35,0.45],1,3)];
g_PS2PNCR=repmat([0.25,0.3,0.35,0.45],1,3);
g_PS2PNCR=[g_PS2PNCR,0,0.25,0.35,0.45,0.05,0.25,0.35,0.45,0.1,0.25,0.35,0.45];
for ii=1:length(params)
    params(ii).g_PS2PNCL=g_PS2PNCL(ii);
    params(ii).g_PS2PNCR=g_PS2PNCR(ii);
end
tic;
%TPNC_matrix = genPatternAry(params, xs, ys, @genTPNCDev);
TPNC_top = genPatternAry(params(1:end/2), xs, ys(1), @genTPNCDev);
TPNC_bot = genPatternAry(params(end/2+1:end), xs, ys(2), @genTPNCDev);
TPNC_top.rotate(180);
TPNC_top.translate([0;381.31]);
toc;

TPNC_matrix=Group(0,0,{TPNC_top,TPNC_bot});
TPNC_matrix.translate([0;-490.22])
%%
fname = 'TPNC03_20190324';
TPNC_matrix.todxf(fname);


