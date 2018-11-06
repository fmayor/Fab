%% generate LNNBZZ
% 42 devices in total
params = genParamStructs('d_nb2refl', [1.75, 1.85],...
    'scale_NB', [1.02:0.005:1.12],...
    'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[3.5],'tether_w_tether', [1],'w2',0.05,...
    'isGenBB',false);
params_ZZ = genParamStructs('g_metal',[0.12],'g_zigzag',0.3,'w_zigzag',[0.3, 0.45],...
    'zigzag_tether', [0.05, 0.1], 'L2', [10,15,20,15]);
params_ZZ = [params_ZZ params_ZZ params_ZZ];
params = copyfield(params, params_ZZ);
% reorder params randomly, generated once and then saved, so that it's
% repeatable
% rnds = rand(1,length(params));
% [~,inds] = sort(rnds);
load('randind_181017.mat');
params = params(inds);
% gen BB once every four devices shield parameters
for ii = 1:length(params)
    if mod(ii,4) == 0
        params(ii).isGenBB = true;
    end
end
dx = 62.5;
nx = length(params);
xs = linspace(-dx*(nx-1)/2,dx*(nx-1)/2, nx);
ys = [0];
x_tot = dx * nx;
fprintf('Generating LNNBZZs...\n');
tic;
LNNBZZ_Row = genPatternAry(params, xs, ys, @test_LNNB_Zigzag);
toc;
LNNBZZ_Row.translate(0, 1700);
LNNBZZ_Row.todxf('tmp');
% devMap.addelement(LNNBZZ_Row);
% LNNBEC_Row_bot = LNNBEC_Row.copy(); LNNBEC_Row_bot.rotate(180);
% devMap.addelement(LNNBEC_Row_bot);