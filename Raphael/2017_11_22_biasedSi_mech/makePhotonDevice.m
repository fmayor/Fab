function photon_dev = makePhotonDevice(ii, chip_tag, hide_GC_holes)

addpath(genpath('/home/tim/github_local/CAD'))
addpath(genpath('/media/user_data/jwitmer/CAD_working_folder/2017_09_26_QuarterWaveQEOM_v1/'))

%% BOOLEANS
B.generate_CAD = true;
B.include_label = false;
B.include_align_marks = false;  %broken
B.invert_tags = false;
B.hide_GC_holes = hide_GC_holes;
B.heidelberg_si_protect = true;

if B.generate_CAD
    P.bragg_point_res = 0.010; % use fine resolution for CAD only, not for plotting
else
    P.bragg_point_res = 1.00;
end

%% LAYERS
% L for Layers
L.SI_LAYER = 'SiLayer';
L.NEG_LAYER = 'NegLayer';
L.HOLE_LAYER = 'HoleLayer'; % Si to be removed from positive regions.
L.INNER_POS_LAYER = 'InnerPosLayer'; % for ring inside
L.HB_SILICON_PROTECT_LAYER = 'HBSiProtectLayer';
L.layer_names = {'SiLayer', 'NegLayer', 'HoleLayer', 'InnerPosLayer','HBSiProtectLayer'};

P.write_field_x = 60;
P.write_field_y = 60;
P.trench = 3.0; % trench around normal waveguides
P.taper_length = 310; % for grating couplers
P.total_wg_length = 520; % length between end of GC tapers (default: 400)
P.tag_offset = [-10;-45];
P.ring_gc_horz_spacing = 46;
P.wg_horz_spacing = P.write_field_x;  %space waveguides by the write field
P.N_points = 3; %isnt used
P.si_protect_width = 1.5;
device.cavity_length = 420;

master_group = Group([0,0],[1;0],[0;1],{});
origin = [0; 0];
curr_x = origin(1);
curr_y = origin(2);

if chip_tag == 'A' || chip_tag == 'B'
  % grating coupler parameters
  P.scaling_pitch = 0.9663;
  P.GC_hole_corr = false;
  % reflector parameters
  top_N_mirror = [4 5 6 9 0];
  device.bot_N_mirror = 500;
  device.bot_N_ramp_in = 10;
  device.bot_N_ramp_out = device.bot_N_ramp_in;
  device.top_N_ramp_in = 10;
  device.top_N_ramp_out = device.top_N_ramp_in;
  device.ampl = 0.110;
  device.period = 0.420;
  device.wg_width = 0.328; % 290 nm
  device.width_corr = false;
  device.ramp_type = 'lin';
  device.type = 'sine';
elseif chip_tag == 'C'
  %dont populate with si devices
  return;
  % grating coupler parameters
  P.scaling_pitch = 1.0;
  P.GC_hole_corr = true;
  % reflector parameters
  top_N_mirror = [6 6 5 3 0];
  device.bot_N_mirror = 500;
  device.bot_N_ramp_in = 10;
  device.bot_N_ramp_out = device.bot_N_ramp_in;
  device.top_N_ramp_in = 10;
  device.top_N_ramp_out = device.top_N_ramp_in;
  device.ampl = 0.100;
  device.period = 0.436;
  device.wg_width = 0.372; % 290 nm
  device.width_corr = true;
  device.duty_cycle = 0.70;
  device.ramp_type = 'lin';
  device.type = 'rect';
end
device.top_N_mirror = top_N_mirror(ii);

GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes, B.heidelberg_si_protect, L.HB_SILICON_PROTECT_LAYER);

photon_dev = makeReflFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B);

end
%make upper left mark to frame the pattern. This helps guide BEAMER's auto-extend feature.
%mk = Rect([-1.5*P.write_field_x + 23; 8.5*P.write_field_y], 0.1, 0.1, L.NEG_LAYER,'base','corner');
%master_group.addelement(mk);


function strip_wg = makeFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  reflector_length = (device.N_ramp_in + device.N_mirror + device.N_ramp_out) * device.period;

  if B.heidelberg_si_protect == true
      total_length = 2*reflector_length + device.cavity_length;
      width = 2*P.trench + device.wg_width - 2*P.si_protect_width;
      elements{end+1} = Rect([0; 0], width, P.total_wg_length, L.HB_SILICON_PROTECT_LAYER);
  end

  %make bottom grating coupler
  elements{end+1} = GC_fn();
  elements{end}.translate([0; 0 - P.total_wg_length/2]);
  elements{end}.rotate([0; 0 - P.total_wg_length/2], -90);

  %make bottom strip waveguide
  strip_length = 0.5 * (P.total_wg_length - device.cavity_length - 2*reflector_length);
  elements{end+1} = makeStripWaveguide(device.wg_width, strip_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([0; -P.total_wg_length/2 + strip_length/2]);

  %make bottom fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.N_ramp_in, device.N_mirror, device.N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.N_ramp_in, device.N_mirror, device.N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);

  %make strip waveguide cavity (at origin)
  elements{end+1} = makeStripWaveguide(device.wg_width, device.cavity_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  %make top fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.N_ramp_in, device.N_mirror, device.N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.N_ramp_in, device.N_mirror, device.N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top strip waveguide
  strip_length = 0.5 * (P.total_wg_length - device.cavity_length - 2*reflector_length);
  elements{end+1} = makeStripWaveguide(device.wg_width, strip_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([0; -P.total_wg_length/2 + strip_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top grating coupler
  elements{end+1} = GC_fn();
  elements{end}.translate([0; 0 - P.total_wg_length/2]);
  elements{end}.rotate([0; 0 - P.total_wg_length/2], -90);
  elements{end}.mirror([0; 0], [1;0]);

  if B.include_label
      tag_pos = [0; 0  - P.total_wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER);
  end
  % make a Group out of the elements
  strip_wg = Group([x0,y0],[1;0],[0;1],elements);
end

function strip_wg = makeAsymFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  bot_reflector_length = (device.bot_N_ramp_in + device.bot_N_mirror + device.bot_N_ramp_out) * device.period;

  top_reflector_length = (device.top_N_ramp_in + device.top_N_mirror + device.top_N_ramp_out) * device.period;

  %make bottom grating coupler
  elements{end+1} = GC_fn();
  elements{end}.translate([0; 0 - P.total_wg_length/2]);
  elements{end}.rotate([0; 0 - P.total_wg_length/2], -90);

  %make bottom strip waveguide
  strip_length = 0.5 * (P.total_wg_length - device.cavity_length - 2*bot_reflector_length);
  elements{end+1} = makeStripWaveguide(device.wg_width, strip_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([0; -P.total_wg_length/2 + strip_length/2]);

  %make bottom fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.bot_N_ramp_in, device.bot_N_mirror, device.bot_N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.bot_N_ramp_in, device.bot_N_mirror, device.bot_N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);

  %make strip waveguide cavity (at origin)
  elements{end+1} = makeStripWaveguide(device.wg_width, device.cavity_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  %make top fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.top_N_ramp_in, device.top_N_mirror, device.top_N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.top_N_ramp_in, device.top_N_mirror, device.top_N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top strip waveguide
  strip_length = 0.5 * (P.total_wg_length - device.cavity_length - 2*top_reflector_length);
  elements{end+1} = makeStripWaveguide(device.wg_width, strip_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([0; -P.total_wg_length/2 + strip_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top grating coupler
  elements{end+1} = GC_fn();
  elements{end}.translate([0; 0 - P.total_wg_length/2]);
  elements{end}.rotate([0; 0 - P.total_wg_length/2], -90);
  elements{end}.mirror([0; 0], [1;0]);

  if B.include_label
      tag_pos = [0; 0  - P.total_wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER);
  end
  % make a Group out of the elements
  strip_wg = Group([x0,y0],[1;0],[0;1],elements);
end

function strip_wg = makeReflFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B)

  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  bot_reflector_length = (device.bot_N_ramp_in + device.bot_N_mirror + device.bot_N_ramp_out) * device.period;

  top_reflector_length = (device.top_N_ramp_in + device.top_N_mirror + device.top_N_ramp_out) * device.period;

  if B.heidelberg_si_protect == true
      total_length = bot_reflector_length + top_reflector_length + device.cavity_length;
      width = 2*P.trench + device.wg_width - 2*P.si_protect_width;
      elements{end+1} = Rect([0; 0], width, P.total_wg_length, L.HB_SILICON_PROTECT_LAYER);
      elements{end+1} = Rect([0; -device.cavity_length/2 - bot_reflector_length/2], width, bot_reflector_length, L.HB_SILICON_PROTECT_LAYER);
  end

  %make bottom fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.bot_N_ramp_in, device.bot_N_mirror, device.bot_N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.bot_N_ramp_in, device.bot_N_mirror, device.bot_N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);

  %make strip waveguide cavity (at origin)
  elements{end+1} = makeStripWaveguide(device.wg_width, device.cavity_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  %make top fishbone reflector
  if device.type == 'sine'
    elements{end+1} = makeFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.top_N_ramp_in, device.top_N_mirror, device.top_N_ramp_out, device.ramp_type, P.bragg_point_res, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  elseif device.type == 'rect'
    elements{end+1} = makeRectFishboneReflector(device.wg_width, P.trench, device.ampl, device.period, device.duty_cycle, device.top_N_ramp_in, device.top_N_mirror, device.top_N_ramp_out, L.SI_LAYER, L.NEG_LAYER, device.width_corr);
  end
  elements{end}.mirror([0;0], [1;0]);
  elements{end}.translate([0;-device.cavity_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top strip waveguide
  strip_length = 0.5 * (P.total_wg_length - device.cavity_length - 2*top_reflector_length);
  elements{end+1} = makeStripWaveguide(device.wg_width, strip_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([0; -P.total_wg_length/2 + strip_length/2]);
  elements{end}.rotate([0; 0], 180);

  %make top grating coupler
  elements{end+1} = GC_fn();
  elements{end}.translate([0; 0 - P.total_wg_length/2]);
  elements{end}.rotate([0; 0 - P.total_wg_length/2], -90);
  elements{end}.mirror([0; 0], [1;0]);

  if B.include_label
      tag_pos = [0; 0  - P.total_wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER);
  end
  % make a Group out of the elements
  strip_wg = Group([x0,y0],[1;0],[0;1],elements);
end

function wg2_ring = makeDoubleWaveguideRingDevice( curr_x, curr_y, GC_fn1, GC_fn2, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

  if B.heidelberg_si_protect == true
      r = device.ring_radius + device.wg_width/2 + P.trench - P.si_protect_width;
      elements{end+1} = Circ([curr_x + device.wg_width + device.ring_cpl_spacing + device.ring_radius; curr_y], r, L.HB_SILICON_PROTECT_LAYER);
  end

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

function thru_wg = makeStripWaveguideDevice( curr_x, curr_y, GC_fn, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

  if B.heidelberg_si_protect == true
      width = 2*P.trench + device.wg_width - 2*P.si_protect_width;
      elements{end+1} = Rect([0; 0], width, P.total_wg_length, L.HB_SILICON_PROTECT_LAYER);
  end

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - device.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - device.wg_length/2], -90);

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - device.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - device.wg_length/2], -90);
  elements{end}.mirror([curr_x;curr_y], [1;0]);

  elements{end+1} = makeStripWaveguide(device.wg_width, device.wg_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  if B.include_label
      tag_pos = [curr_x; curr_y  - device.wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER);
  end
  % make a Group out of the elements
  thru_wg = Group([x0,y0],[1;0],[0;1],elements);
end
