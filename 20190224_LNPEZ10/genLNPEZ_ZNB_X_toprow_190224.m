function [LNPEZ_ZNB_X_Row, L_ZNBX] = genLNPEZ_ZNB_X_toprow_190224()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz
%% params for ZNBZ
% paramsZNB1 = genParamStructs('d_nb2refl', [2.2,2.3],...
%     'scale_NB', [0.89:0.005:1.0],...
%     'w_end', 0.26, 'tether_w2', [2.4],...
%     'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
%     'isGenBB',false, 'isGenPS',false,'theta',-90);
% 
% paramsZNB2 = genParamStructs('L',[5,10,14],'L2',0,'g',0.3,'w',0.45,'g_metal',0.15);
% paramsZNB2 = repmat([paramsZNB2 circshift(paramsZNB2,1) circshift(paramsZNB2,2)],1,15);
% paramsZNB = copyfield(paramsZNB1, paramsZNB2);
% 
% gaps = [0.15,0.16,0.17,0.18,0.19];
% gaps = num2cell(repmat(gaps, 1, 40));
% [paramsZNB.gap] = gaps{1:length(paramsZNB)};
%% params for LNNBZZ (tunable kappa_e)
% % 42 devices in total
% paramsNBZZ = genParamStructs('scale_NB', [1.0:0.0033:1.08],...
%     'w_end', 0.26, 'tether_w2', [2.4],...
%     'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
%     'isGenBB',false);
% params_ZZ = genParamStructs('d_nb2refl', [1.8,1.85,1.9,1.95,2.0,2.05],...
%     'd_nb2refl2', [1.75, 1.85, 1.75, 1.85, 1.85, 1.75,1.85, 1.75],'g_metal',[0.15],...
%     'g_zigzag',0.3,'w_zigzag',[0.45],...
%     'zigzag_tether', [0.05, 0.1, 0.05, 0.1, 0.1, 0.05, 0.1, 0.05], 'L2', [10,15,10,15,15,10,15,10]);
% params_ZZ = [params_ZZ circshift(params_ZZ,1) params_ZZ circshift(params_ZZ,1)];
% paramsNBZZ = copyfield(paramsNBZZ, params_ZZ);
% % reorder params randomly, generated once and then saved, so that it's
% % repeatable
% % rnds = rand(1,length(paramsNBZZ));
% % [~,inds] = sort(rnds);
% % load('randind_181017.mat');
% % paramsNBZZ = paramsNBZZ(inds);
% % gen BB once every four devices shield parameters
% for ii = 1:length(paramsNBZZ)
%     if mod(ii,4) == 0
%         paramsNBZZ(ii).isGenBB = true;
%     end
% end
%% params for PEZ
% paramsPEZ1 = genParamStructs('w_mid', [0.5,0.55,0.6,0.65],...
%     'w_end', 0.26,...
%      'g_mid', [0.3],'zz_w',[0.45],'g_metal',[0.15],...
%      'offset',[0,0.2,0.3,0.35,0.4,0.5,0.6]);
%  
% paramsPEZ2 = genParamStructs('L', [10,15,20],'L2',[5,10],'g', [0.3,0.6,0.3,0.6,0.6,0.3,0.6,0.3]);
% paramsPEZ2 = [paramsPEZ2 paramsPEZ2 paramsPEZ2];
% paramsPEZ = copyfield(paramsPEZ1, paramsPEZ2);

%% params for ZNB X
paramsZNBX1 = genParamStructs('scale_NB', [0.89:0.0054:1.08],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false,'theta',-90);

paramsZNBX2 = genParamStructs('gap',[0.23,0.28,0.32],'L',[8,12,17,22],'d_nb2refl', [2,2.1,2.2],'L2',0,'g',0.3,...
    'w',0.45,'g_metal',0.15);
paramsZNBX2 = repmat([paramsZNBX2 circshift(paramsZNBX2,1) circshift(paramsZNBX2,2)],1,15);
paramsZNBX = copyfield(paramsZNBX1, paramsZNBX2);


% gaps = [0.15,0.17,0.2,0.22,0.25,0.27,0.29,0.32,0.35,0.38,0.41];
% gaps = num2cell(repmat(gaps, 1, 40));
% [paramsZNBX.gap] = gaps{1:length(paramsZNBX)};


%% params for rotated ZNB X (bad design on purpose)
paramsZNBX1 = genParamStructs('scale_NB', [0.89:0.0215:1.08],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false,'theta',0);

paramsZNBX2 = genParamStructs('gap',[0.28],'L',[5,10,20],'d_nb2refl', [2,2.1,2.2],'L2',0,'g',0.3,...
    'w',0.45,'g_metal',0.15);
paramsZNBX2 = repmat([paramsZNBX2 circshift(paramsZNBX2,1) circshift(paramsZNBX2,2)],1,15);
paramsZNBX0 = copyfield(paramsZNBX1, paramsZNBX2);

%%
dx0=2*62.5;
nx0 = length(paramsZNBX0);
xs0 = linspace(-dx0*(nx0-1)/2,dx0*(nx0-1)/2, nx0);
ys0 = [0];
L_ZNBX0 = dx0*nx0;

dx2 = 2*62.5;
nx2 = length(paramsZNBX);
xs2 = linspace(-dx2*(nx2-1)/2,dx2*(nx2-1)/2, nx2);
ys2 = [0];
L_ZNBX = dx2*nx2;

L=L_ZNBX0+L_ZNBX;
fprintf('Generating LNZNBs...\n');
tic;
LNZNBX_Row = genPatternAry(paramsZNBX, xs2, ys2, @genLNZNB_DS_EC_ZZ_X);
LNZNBX_Row.translate([-L/2+L_ZNBX/2;0]);

LNZNBX0_Row = genPatternAry(paramsZNBX0, xs0, ys0, @genLNZNB_DS_EC_ZZ_X);
%LNZNB_X_Row.translate([L_ZNBX/4+L_PEZ/4;0]);
LNZNBX0_Row.translate([-L/2+L_ZNBX+L_ZNBX0/2+55;0]);
toc;
%% fix problems manually
% %feed
% feed_link=Rect(-1219,-1062,15,13,'base','center');
% feed_link.layer = 'metal_BB';
%%
LNPEZ_ZNB_X_Row=gpack.Group(0,0,{LNZNBX_Row,LNZNBX0_Row});


end