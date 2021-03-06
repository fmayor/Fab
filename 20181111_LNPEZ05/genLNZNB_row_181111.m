function [LNZNB_Row, L_NB] = genLNZNB_row_181111()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz

params1 = genParamStructs('d_nb2refl', [2.2,2.3],...
    'scale_NB', [0.89:0.006:1.0],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false,'theta',-90);

params2 = genParamStructs('L',[5,10,14],'L2',0,'g',0.3,'w',0.45,'g_metal',0.15);
params2 = repmat([params2 circshift(params2,1) circshift(params2,2)],1,15);
params = copyfield(params1, params2);


gaps = [0.15,0.16,0.17,0.18,0.19];
gaps = num2cell(repmat(gaps, 1, 40));
[params.gap] = gaps{1:length(params)};

dx = 2*62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
L_NB = dx*nx;
fprintf('Generating LNZNBs...\n');
tic;
LNZNB_Row = genPatternAry(params, xs, ys, @genLNZNB_DS_EC_ZZ_X);
toc;

end