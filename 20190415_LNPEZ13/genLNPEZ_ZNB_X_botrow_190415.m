function [LNPEZ_bot_Row, L_LNR] = genLNPEZ_ZNB_X_botrow_190415()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz

%% params for LNR
paramsLNR = genParamStructs('scale_NB', [0.9:0.02:1.1],...
    'ring_d', [0.2:0.1:0.7],'w_end',0.26,'ring_radius',20, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false,'theta',0,'d_nb2refl', [0.1],'d_nb2refl2', [0.1]);
refl_offset = [20,30,20,30,20,30,30,20,30,20,30,20];
refl_offset = num2cell(repmat(refl_offset, 1, 40));
[paramsLNR.refl_offset] = refl_offset{1:length(paramsLNR)};
%% params for ZNB X
% paramsZNBX1 = genParamStructs('scale_NB', [0.93:0.0024:1.01],...
%     'w_end', 0.26, 'tether_w2', [2.4],...
%     'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
%     'isGenBB',false, 'isGenPS',false,'theta',-90);
% 
% paramsZNBX2 = genParamStructs('gap',[0.28,0.35],'L',[2,6,8,10,14,16,19,20],'d_nb2refl', [2.1,2.2,2.3],'L2',0,'g',0.3,...
%     'w',0.45,'g_metal',0.15);
% paramsZNBX2 = repmat([paramsZNBX2 circshift(paramsZNBX2,1) circshift(paramsZNBX2,2)],1,15);
% paramsZNBX = copyfield(paramsZNBX1, paramsZNBX2);


% gaps = [0.15,0.17,0.2,0.22,0.25,0.27,0.29,0.32,0.35,0.38,0.41];
% gaps = num2cell(repmat(gaps, 1, 40));
% [paramsZNBX.gap] = gaps{1:length(paramsZNBX)};
%%
dx = 90;%2*62.5;
nx = length(paramsLNR);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
L_LNR = dx*nx;

L=L_LNR;
fprintf('Generating Rings...\n');
tic;
% LNNBZZ_Row = genPatternAry(paramsNBZZ, xs0, ys0, @test_LNNB_Zigzag);
% LNNBZZ_Row.translate([-L/2+L_NBZZ/2;0]);
% 
% LNPEZ_Row = genPatternAry(paramsPEZ, xs, ys, @genTvsDx_withZZ);
% %LNPEZ_Row.translate([-L_ZNBX/4-L_PEZ/4;0]);
% LNPEZ_Row.translate([-L/2+L_NBZZ+L_PEZ/2;0]);

LNR_Row = genPatternAry(paramsLNR, xs, ys, @genRingCoupler);
%LNZNB_X_Row.translate([L_ZNBX/4+L_PEZ/4;0]);
%LNR_Row.translate([-L/2-L_LNR/2;0]);
toc;
%% fix problems manually
%feed
% feed_link=Rect(-1219,-1062,15,13,'base','center');
% feed_link.layer = 'metal_BB';
%%
LNPEZ_bot_Row=gpack.Group(0,0,{LNR_Row});


end