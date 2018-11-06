% test LNZNB

%%

global config;
% config.isfast = true;
config.isfast = false;
config.mkrPixelSize = 0.6;

param = struct('scale_NB', 1,...
    'd_nb2refl', 1.8,'w_end', 0.26, 'tether_w2', [2.4],...
    'tether_L',[4],'tether_w_tether', [1],'isGenPS',true,'w2',0.05);
param.ind = 0;
dev4SEM = genLNZNB_DS_EC(param);
dev4SEM.todxf('test_LNZNB');
