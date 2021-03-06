function [LNPEZ_ZNB_Z_Row, L_ZNB] = genLNPEZ_ZNB_Z_row_181111()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz

%params for ZNB
paramsZNB1 = genParamStructs('d_nb2refl', [2.2,2.3],...
    'scale_NB', [0.89:0.005:1.0],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false,'theta',-90);

paramsZNB2 = genParamStructs('L',[5,10,14],'L2',0,'g',0.3,'w',0.45,'g_metal',0.15);
paramsZNB2 = repmat([paramsZNB2 circshift(paramsZNB2,1) circshift(paramsZNB2,2)],1,15);
paramsZNB = copyfield(paramsZNB1, paramsZNB2);

gaps = [0.15,0.16,0.17,0.18,0.19];
gaps = num2cell(repmat(gaps, 1, 40));
[paramsZNB.gap] = gaps{1:length(paramsZNB)};

%params for PEZ
paramsPEZ1 = genParamStructs('w_mid', [0.5,0.55,0.6,0.65],...
    'w_end', 0.26,...
     'g_mid', [0.15],'zz_w',[0.45],'g_metal',[0.15],...
     'offset',[0,0.2,0.3,0.35,0.4,0.5,0.6]);
 
paramsPEZ2 = genParamStructs('L', [10,15,20],'L2',[5,10,15],'g', [0.3,0.6]);
paramsPEZ2 = [paramsPEZ2 paramsPEZ2];
paramsPEZ = copyfield(paramsPEZ1, paramsPEZ2);

dx = 62.5;
nx = length(paramsZNB);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
L_ZNB = dx*nx;
dx2 = 62.5;
nx2 = length(paramsPEZ);
xs2 = linspace(-dx2*(nx2-1)/2,dx2*(nx2-1)/2, nx2);
ys2 = [0];
L_PEZ = dx2*nx2;
fprintf('Generating LNZNBs...\n');
tic;
LNZNB_Z_Row = genPatternAry(paramsZNB, xs, ys, @genLNZNB_DS_EC_ZZ_Z);
LNZNB_Z_Row.translate([-L_ZNB/4-L_PEZ/4;0]);
LNPEZ_Row = genPatternAry(paramsPEZ, xs2, ys2, @genTvsDx_withZZ);
LNPEZ_Row.translate([L_ZNB/4+L_PEZ/4;0]);
toc;
LNPEZ_ZNB_Z_Row=gpack.Group(0,0,{LNZNB_Z_Row, LNPEZ_Row});

end