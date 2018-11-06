function [rng_group, dev_counter] = ring_sweep(P, L, B, dev_counter)

addpath(genpath('/home/tim/github_local/CAD'))

P.scaling_pitch = 0.9663;
P.GC_hole_corr = false;

rng_group = Group([0,0],[1;0],[0;1],{});
origin = [0; 0];
curr_x = origin(1);
curr_y = origin(2);

dev.wg_width = 0.300;
ring_cpl_spacing = repmat([0.400:0.100:1.0], 1, 2);
ring_gc_horz_spacing = [46*ones(1,7) 92*ones(1,7)];

for ii = 1:length(ring_cpl_spacing)
  if dev_counter == 79
    curr_x = curr_x + 37;
  end
  dev.wg_length = 150;
  dev.ring_cpl_spacing = ring_cpl_spacing(ii);
  P.ring_gc_horz_spacing = ring_gc_horz_spacing(ii);
  dev.ring_radius = P.ring_gc_horz_spacing/2 - dev.wg_width - dev.ring_cpl_spacing;

  dev.tag_text1 = num2str(dev_counter);
  dev.tag_text2 = num2str(dev_counter + 1);

  GC_fn = @() makeTwoStepGrating(dev.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);

  ringres_dev = makeDoubleWaveguideRingDevice(curr_x, curr_y, GC_fn, GC_fn, dev, P, L, B);
  rng_group.addelement(ringres_dev);
  dev_counter = dev_counter + 2;
  curr_x = curr_x + 2*P.write_field_x;
end

end

function wg2_ring = makeDoubleWaveguideRingDevice( curr_x, curr_y, GC_fn1, GC_fn2, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

  %make left waveguide
  device.tag_text = device.tag_text1;
  wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn1, device, P, L, B);
  elements{end+1} = wg_dev;

  %make ring resonator
  elements{end+1} = RingRes([curr_x + device.wg_width + device.ring_cpl_spacing + device.ring_radius; curr_y], device.ring_radius, device.wg_width, P.trench, P.trench, P.trench, P.N_points, L.SI_LAYER, L.NEG_LAYER, L.HOLE_LAYER, L.INNER_POS_LAYER);

  if B.include_label
      tag_pos = [curr_x; curr_y  - device.wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text1, B, L.NEG_LAYER);
  end

  %make right waveguide
  device.tag_text = device.tag_text2;
  curr_x = curr_x + P.ring_gc_horz_spacing;
  wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn2, device, P, L, B);
  elements{end+1} = wg_dev;

  % make a Group out of the elements
  wg2_ring = Group([x0,y0],[1;0],[0;1],elements);
end

function thru_wg = makeStripWaveguideDevice( curr_x, curr_y, GC_fn, dev, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

  if B.heidelberg_si_protect == true
      width = 2*P.trench + dev.wg_width - 2*P.si_protect_width;
      elements{end+1} = Rect([0; 0], width, P.total_wg_length, L.HB_SILICON_PROTECT_LAYER);
  end

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - dev.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - dev.wg_length/2], -90);

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - dev.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - dev.wg_length/2], -90);
  elements{end}.mirror([curr_x;curr_y], [1;0]);

  elements{end+1} = makeStripWaveguide(dev.wg_width, dev.wg_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  if B.include_label
%      tag_pos = [curr_x; curr_y  - dev.wg_length/2 - P.taper_length ] + P.tag_offset ;
  tag_pos = [curr_x; curr_y  - dev.wg_length/2 - P.taper_length] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, dev.tag_text, B, L.NEG_LAYER);
  end
  % make a Group out of the elements
  thru_wg = Group([x0,y0],[1;0],[0;1],elements);
end
