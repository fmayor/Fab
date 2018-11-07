function [LNZNB_Row, L_NB] = genLNZNB_row_181106()


% center the scaling to be 96% since the fundamental mode is at 184.5 THz

params = genParamStructs('d_nb2refl', [1.8, 2, 2.2],...
    'scale_NB', [0.88:0.006:1.03],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false, 'isGenPS',false );
gaps = [0.14, 0.15, 0.16, 0.17, 0.18, 0.19];
gaps = num2cell(repmat(gaps, 1, 15));
[params.gap] = gaps{1:length(params)};

dx = 62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
L_NB = dx*nx;
fprintf('Generating LNZNBs...\n');
tic;
LNZNB_Row = genPatternAry(params, xs, ys, @genLNZNB_DS_EC);
toc;

end