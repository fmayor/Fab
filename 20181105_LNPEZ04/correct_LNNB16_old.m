function dev = correct_LNNB16(dev)
% measured actual size divided by designed size from SEM
c_mir.wmax = 0.9816;
c_mir.wmin = 1.0051;
c_mir.hy = 0.9318;
c_mir.hx = 0.9403;

c_def.wmax = 0.9828;
c_def.wmin = 0.9957;
c_def.hy = 0.9287;
c_def.hx = 0.8984;

c_w_end = 0.9217;

% apply correction measured from LNNB16, 20180828
dev.w_end = dev.w_end / c_w_end;
dev.P_mirror.hx = dev.P_mirror.hx / c_mir.hx;
dev.P_mirror.hy = dev.P_mirror.hy / c_mir.hy;

dev.P_defect.hx = dev.P_defect.hx / c_def.hx;
dev.P_defect.hy = dev.P_defect.hy / c_def.hy;

dev.P_mirror = correct_W_amp(dev.P_mirror, c_mir.wmax, c_mir.wmin);
dev.P_defect = correct_W_amp(dev.P_defect, c_def.wmax, c_def.wmin);

% correction measured on 180925, still from LNNB16
try
    dev.P_1DPS.w1x = dev.P_1DPS.w1x/0.9699;
    dev.P_1DPS.w1y = dev.P_1DPS.w1y/0.9913;
catch err
    
end
end

function P = correct_W_amp(P, corr_wmax, corr_wmin)

w_max = P.w + 2 * P.amp;
w_min = P.w - 2 * P.amp;
w_max_c = w_max/corr_wmax;
w_min_c = w_min/corr_wmin;
w_c = (w_max_c + w_min_c)/2;
amp_c = (w_max_c - w_min_c)/4;
P.w = w_c;
P.amp = amp_c;

end

