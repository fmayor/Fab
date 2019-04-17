% generate TPNC devices
% WTJ,FMM, 20190415


%%
% I have 5000 by 2200 um, could fit in ~ 12 by 2 TPNC devices
dx = 400;
dy = 600;

nx = 12;
ny = 2;
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = linspace(-dy*(ny-1)/2,dy*(ny-1)/2, ny);
%meshgrid is done within genPatternAry
%% top row
tic;
params = genParamStructs('w_def', [1.12],'L_Bender',[20,25,30],...
'g_PS2PNCL',[0,0,0,0],'flag_def',false); 
%make different gaps for left and right tips
g_PS2PNCL=[0.1,0.2,0.35,0.1,0.15,0.25,0.3,0.3,0.1,0.25,0.35,0.15];
g_PS2PNCR=[0.25,0.25,0.1,0.25,0.3,0.25,0.15,0.1,0.3,0.25,0.35,0.3];
for ii=1:length(params)
    params(ii).g_PS2PNCL=g_PS2PNCL(ii);
    params(ii).g_PS2PNCR=g_PS2PNCR(ii);
end
TPNC_top = genPatternAry(params(1:end), xs, ys(1), @genTPNCDev);
TPNC_top.rotate(180);
TPNC_top.translate([0;381.31]);
%% bot row
params = genParamStructs('w_def', [1.12],'L_Bender',[20,25,30],...
'g_PS2PNCL',[0,0,0,0],'flag_def',true); 
%make different gaps for left and right tips
g_PS2PNCL=[0.1,0.15,0.25,0.3,0.1,0.2,0,0.25,0.15,0.3,0.2,0.1];
g_PS2PNCR=g_PS2PNCL;
for ii=1:length(params)
    params(ii).g_PS2PNCL=g_PS2PNCL(ii);
    params(ii).g_PS2PNCR=g_PS2PNCR(ii);
end
TPNC_bot = genPatternAry(params(1:end), xs, ys(2), @genTPNCDev);

%%

TPNC_matrix=Group(0,0,{TPNC_top,TPNC_bot});
TPNC_matrix.translate([0;-490.22])
%%
fname = 'TPNC04_20190415';
TPNC_matrix.todxf(fname);


