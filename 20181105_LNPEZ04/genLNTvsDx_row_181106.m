function [LNTvsDx_row, L_NBZZ] = genLNTvsDx_row_181106()



params_passive = genParamStructs('w_mid', [0.3,0.4,0.5,0.6],...
    'w_end', 0.26,'g_mid', 0.15,...
    'offset',[0:0.1:0.4,0.45,0.5,0.55,0.6:0.1:1,1.5,2.0]);
     
 
params_active1 = genParamStructs('w_mid', [0.4,0.5,0.6],...
    'w_end', 0.26,...
     'g_mid', [0.15],'zz_w',[0.45],'g_metal',[0.15],...
     'offset',[0,0.2,0.3,0.4,0.5,0.6,0.8]);
 
params_active2 = genParamStructs('L', [10,15,20],'L2',[10,15],'g', [0.3,0.6]);
params_active2 = [params_active2 params_active2];
params_active = copyfield(params_active1, params_active2);


dx = 62.5;
nx = length(params_passive)+length(params_active);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
xs2=xs(1:4:4*length(params_active));
ys2=[0];
xs1=setdiff(xs,xs2);
ys1=[0];

L_NBZZ = dx*nx;
fprintf('Generating LNTvsDx...\n');
tic;
LNTvsDx_passive = genPatternAry(params_passive, xs1, ys1, @genTvsDx_noZZ_three_devs);
LNTvsDx_active = genPatternAry(params_active, xs2, ys2, @genTvsDx_withZZ);
toc;
LNTvsDx_row = gpack.Group(0,0,{});
LNTvsDx_row.addelement(LNTvsDx_passive);
LNTvsDx_row.addelement(LNTvsDx_active);
end