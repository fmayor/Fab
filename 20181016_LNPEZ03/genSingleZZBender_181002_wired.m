function [dev,port_L1,port_R1,port_L2,port_R2] = genSingleZZBender_181002_wired(param)
% two zzb connected by a short piece of WG with joints
% WTJ, 20180628
% 
% Modified for PEZZS
% WTJ & FM, 181002

param_singleZZB.w = param.w;
param_singleZZB.g = param.g;
param_singleZZB.L2 = param.L2;
param_singleZZB.L = param.L;
param_singleZZB.g_metal=param.g_metal;
param_singleZZB.layer_sharp = param.layer_sharp;
param_singleZZB.layer_metalBB = param.layer_metalBB;
param_singleZZB.layer_metalSB = param.layer_metalSB;
param_singleZZB.layer_uct = param.layer_uct;
w_joint = param.w_joint;
l_joint = 0.25;
if isfield(param, 'l_joint')
    l_joint = param.l_joint;
end
L_WG = 2; % WG length between two zzb
w_WG = 0.8;

%if param.flag_signFlip==true
    param_singleZZB.flip=false;
    [zzb1,port_L1,port_R1] = genZZBender_longArm_181002_wired(param_singleZZB);
    param_singleZZB.flip=true;
    [zzb2,port_L2,port_R2] = genZZBender_longArm_181002_wired(param_singleZZB);
% else
%     param_singleZZB.flip=false;
%     zzb1 = genZZBender_longArm_181002_wired(param_singleZZB);
%     zzb2 = genZZBender_longArm_181002_wired(param_singleZZB);
% end
zzb2.mirror([0;0],[1;0]);
% top coordinates of zzb1
x_top = mean(zzb1.elements{1}.port2(1,:));
y_top = zzb1.elements{1}.port2(2,1);

zzb1.translate(0,  - y_top);
zzb2.translate(0,  + y_top);
zzb2.rotate(180);
% translate zzb according to zzb's L and the tether length, also traslate
% upward so that y = 0 is aligned to the center of the last arm
d_tetherCenter2ZZEnd = 3.5;  % hard coded, from genPEZZS
dx = d_tetherCenter2ZZEnd + param.L/2 + l_joint + param.w/2;
dy = param.g/2 + param.w/2;
zzb1.translate(dx, dy);
zzb2.translate(-dx, dy);
port_L1=port_L1+[0+dx;-y_top+dy];port_R1=port_R1+[0+dx;-y_top+dy];
port_L2=port_L2+[0+dx;-y_top+dy];port_R2=port_R2+[0+dx;-y_top+dy];
% add joints
r_j1 = Rect(d_tetherCenter2ZZEnd + param.w/2, 0, l_joint * 1.5, w_joint, 'base','center');
r_j1.layer = 'M1_Neg';
r_j2 = r_j1.copy();
r_j2.translate(-2*d_tetherCenter2ZZEnd - param.w, 0);
buff_cure = 3;

dev = gpack.Group(0,0,{zzb1,zzb2,r_j1,r_j2});

% add writefield
% w_F = 60;
% r_RF = Rect(0,0,w_F,w_F,'base','center');
% r_RF.layer = 'M1_Field';
% dev.addelement(r_RF);

% add release window covering WG
w_rls = 13;
r_rls = Rect(0,0,w_rls,w_rls,'base','center');
r_rls.layer = 'M2_undercutMsk';
dev.addelement(r_rls);


% add label
% [ind_pos, ind_neg] = genStringPolygon(sprintf('%d',param.ind),...
%         0.4);
% ind_pos.layer = 'M1_Neg';
% ind_neg.layer = 'M1_Pos_ind';
% ind = gpack.Group(0,0,{ind_pos, ind_neg});
% dx = 0;
% dy = -50;
% ind.rotate(-90);
% ind.translate(dx,dy);
% 
% dev.addelement(ind);

end