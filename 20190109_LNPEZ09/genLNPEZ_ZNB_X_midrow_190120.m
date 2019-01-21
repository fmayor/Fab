function [CL_Row, L_CL] = genLNPEZ_ZNB_X_midrow_190120()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz


%% params for CL
 
paramsCL1 = genParamStructs('L', [5,8,10,13],'L2',[0,1,2],'g', [0.3],'g_metal',[0.15],'w',[0.45]);
paramsCL2 = genParamStructs('L', [3,6,9,15],'L2',[0,1],'g', [0.6],'g_metal',[0.15],'w',[0.45]);
paramsCL=[paramsCL1 paramsCL2];
%%


dx = 250;
nx = length(paramsCL);
xs = [-dx*(nx-1)/2-dx:dx:-3*dx/2,dx/2:dx:dx*(nx-1)/2];
ys = [0];
L_CL = dx*nx;




tic;

CL_Row = genPatternAry(paramsCL, xs, ys, @genCantilever);

toc;
%% fix problems manually
%feed
% feed_link=Rect(-1219,-1062,15,13,'base','center');
% feed_link.layer = 'metal_BB';
%%
%CL_Row=gpack.Group(0,0,{LNPEZ_Row,LNZNB_X_Row,LNNBZZ_Row,feed_link});


end