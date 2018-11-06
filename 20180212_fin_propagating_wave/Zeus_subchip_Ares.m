% Zeus subchip Ares layout
% RVL 2017-11-30
% EOS6 1000pA

function out = Zeus_subchip_Ares(B, chip_tag)

addpath(genpath('C:\Users\rvanlaer\Documents\GitHub\CAD\'))
addpath(genpath('C:\Users\rvanlaer\Documents\GitHub\CAD\device\photonics\RVL_various_220SiPhot\'))

% clear
tic

% dxf file name, without .dxf
filename = 'Zeus_subchip_Ares';
B.heidelberg_si_protect = true;

%% LAYERS
% L for Layers
L.SI_LAYER = 'SiLayer';
L.NEG_LAYER = 'NegLayer';
L.HOLE_LAYER = 'HoleLayer'; % Si to be removed from positive regions.
L.INNER_POS_LAYER = 'InnerPosLayer'; % for ring inside
L.HB_SILICON_PROTECT_LAYER = 'HBSiProtectLayer';% L.METAL_ON_SI_LAYER = 'MetalOnSiLayer'; % Si to be removed from positive regions.
L.HB_METAL_GAP_LAYER = 'HBMetalGapLayer';
L.HB_METAL_POS_LAYER = 'HBMetalPosLayer';
L.EOS6_METAL_LAYER = 'EOS6MetalLayer';
L.DICE_MARK_LAYER = 'DiceMarkLayer';
L.HB_UNDERCUT = 'HBUndercut';

%% Define global parameters
% P for Parameters
P.write_field_x = 60;
P.write_field_y = 60;
P.wg_width = 0.474; % corrected, should come out around 450 nm.
P.bragg_point_res = 0.010; % use fine resolution for CAD only, not for plotting
P.trench = 3; % trench around normal waveguides
P.taper_length = 35; % for grating couplers
P.total_wg_length = 270+290-14; % length between end of GC tapers (default: 400)
P.scaling_pitch = 1.0;  % (-) Scaling of the pitch
P.scaling_holes = 1.0; % (-) Scaling of the grating holes on top of potential correction factor
P.tag_offset = [-10;-40-14];
P.wg_horz_spacing = 480; %P.write_field_x;  %space waveguides by the write field
P.fishbone_wg_width = P.wg_width;
P.ring_gc_horz_spacing = 46;
P.si_protect_width = 1.5; % half the trench is optimal
P.N_points = 3;  %this is for the RingRes function. NEED TO REMOVE FROM FNC!
P.contact_wire_width = 1; % width of wire protruding from slot mode converter
P.contact_wire_length = 80;
P.array_distance = P.write_field_x;
%P.GC_hole_corr = false;
P.electrode_offset = 1.5;

%% Initialize
master_pattern_subchip    =  Group([0,0],[1;0],[0;1],{});

% pattern must be away from [0,0] to avoid errors in AutoCAD, so offset
% origin
origin = [0;0];

curr_x = origin(1) + P.wg_horz_spacing/2;
curr_y = origin(2);

device_counter = 1;
wg_electrode_locations = {};
wg_electrode_total_length = {};
fishbone_total_length = []; % length of cavity, including reflectors, for each fishbone device.

%% Make straight waveguides without electrodes: play with scaling of the holes
N_wgs = 5;
grating_holes_grp = Group([0,0],[1;0],[0;1],{});

%DEVICE PARAMETERS
wg_width = P.wg_width;
% reflector parameters
scaling_holes = repmat([1 0.97 1.03 0.94 1.05], 1, 1);

device = {};

for ii = 1:N_wgs
    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.width_corr = true;
    device.wg_length = P.total_wg_length;
    device.tag_text = num2str(device_counter);
    device.scaling_holes = scaling_holes(ii);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,device.scaling_holes,P.wg_width,P.si_protect_width);
    
    wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    grating_holes_grp.addelement(wg_dev);

    curr_x = curr_x + P.write_field_x;
    device_counter = device_counter + 1;

end
curr_x = curr_x + P.array_distance;
clear device;
master_pattern_subchip.addelement(grating_holes_grp);

%% Make straight waveguides without electrodes: play with focusing length
N_wgs = 4;
grating_focus_grp = Group([0,0],[1;0],[0;1],{});

%DEVICE PARAMETERS
wg_width = P.wg_width;
% reflector parameters
focusing_length = 37.8*repmat([0.9 1.1 0.8 1.2], 1, 1);

device = {};

for ii = 1:N_wgs
    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.width_corr = true;
    device.wg_length = P.total_wg_length;
    device.tag_text = num2str(device_counter);
    device.focusing_length = focusing_length(ii);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,P.wg_width,P.si_protect_width,'focusing_length',device.focusing_length);
    
    wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    grating_focus_grp.addelement(wg_dev);

    curr_x = curr_x + P.write_field_x;
    device_counter = device_counter + 1;

end
curr_x = curr_x + P.array_distance;
clear device;
master_pattern_subchip.addelement(grating_focus_grp);

%% Make straight waveguides without electrodes: play with defocus
N_wgs = 4;
grating_defocus_grp = Group([0,0],[1;0],[0;1],{});

%DEVICE PARAMETERS
wg_width = P.wg_width;
% reflector parameters
defocus = repmat([5 -5 10 -10], 1, 1); % 5 defocus is hero device?!

device = {};

for ii = 1:N_wgs
    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.width_corr = true;
    device.wg_length = P.total_wg_length;
    device.tag_text = num2str(device_counter);
    device.defocus = defocus(ii);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,P.wg_width,P.si_protect_width,'defocus',device.defocus);
    
    wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    grating_defocus_grp.addelement(wg_dev);

    curr_x = curr_x + P.write_field_x;
    device_counter = device_counter + 1;

end
curr_x = curr_x + P.array_distance;
clear device;
master_pattern_subchip.addelement(grating_defocus_grp);

%% Make ring resonator
ring_grp = Group([0,0],[1;0],[0;1],{});

curr_x = curr_x + P.write_field_x - 5;

%Make a standard Si strip ring resonator
device = {};

N_rings = 1;
for ii = 1:N_rings

    device.wg_length = P.total_wg_length;
    device.wg_width = P.wg_width;
    device.ring_cpl_spacing = 0.300;
    device.ring_radius = P.ring_gc_horz_spacing/2 - device.wg_width - device.ring_cpl_spacing ;

    device.tag_text1 = num2str(device_counter);
    device.tag_text2 = num2str(device_counter + 1);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,device.wg_width,P.si_protect_width);
    
    ring_grp.add( makeDoubleWaveguideRingDevice(curr_x, curr_y, GC_fn, GC_fn, device, P, L, B) );
    device_counter = device_counter + 2;
    curr_x = curr_x + 2*P.write_field_x;
end
clear device

ring_grp.translate([-P.ring_gc_horz_spacing/2;0]); % translate ring group back to center ring in write field.
master_pattern_subchip.addelement(ring_grp);

%% Make the waveguide array: play with waveguide width, length, electrode spacing, electrode offset, create a fishbone too

N_wgs = 8;
wg_electrode_grp = Group([0,0],[1;0],[0;1],{});

% want to shift center wavelength from 1540 to 1578, so new period is 434
% nm

%DEVICE PARAMETERS
wg_widths = 1.0533*repmat([0.45 0.5 0.55 0.6 0.45 0.5 0.55 0.45], 1, 1);
wg_length = repmat([546 546 546 546 900 900 900 450], 1, 1);
% reflector parameters
cavity_length = 450;
% reflector parameters
N_ramp_in = 10;
N_mirror = 6; % only the last value is used
N_ramp_out = 10;
ampl = 0.136;
neff = 2.255; % (-) Effective index for 450 nm wide silicon waveguide at 1550 nm
lambda0 = 1.55; % target center wavelength
period = lambda0/(2*neff); % Should be about 344 nm

undercut_length = 15;
undercut_period = 25;
undercut_spacing = 1.25; % (um) Undercut spacing from the waveguide
undercut_widths = wg_widths + 2*undercut_spacing;
number_of_suspensions = floor(wg_length/undercut_period);

device = {};
curr_x = curr_x - P.wg_horz_spacing/2;
for ii = 1:N_wgs

    curr_x = curr_x + P.wg_horz_spacing;

    %DEVICE PARAMETERS
    device.wg_width = wg_widths(ii);
    device.wg_length = wg_length(ii);
    device.tag_text = num2str(device_counter);
    device.width_corr = true;
    
    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,device.wg_width,P.si_protect_width);
    
    if ii < 8
        wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
        for kk = 1:number_of_suspensions(ii)
            suspension = Rect([0; 0], undercut_widths(ii), undercut_length, L.HB_UNDERCUT);
            suspension.translate(0,(kk-1/2)*undercut_period - device.wg_length/2);
            wg_dev.add(suspension);
        end
    end
    %fishbone_dev = makeFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    
    
    if ii == 8
        device.N_ramp_in = N_ramp_in;
        device.N_mirror = N_mirror;
        device.N_ramp_out = N_ramp_out;
        device.ampl = ampl;
        device.period = period;
        device.ramp_type = 'lin';
        device.type = 'rect';
        device.duty_cycle = 0.70;
        device.width_corr = true;
        device.cavity_length = cavity_length;
        wg_dev = makeFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B);
        
        for kk = 1:number_of_suspensions(ii)
            suspension = Rect([0; 0], undercut_widths(ii), undercut_length, L.HB_UNDERCUT);
            suspension.translate(0,(kk-1/2)*undercut_period - device.wg_length/2);
            wg_dev.add(suspension);
        end
    end
    wg_electrode_grp.addelement(wg_dev);

    wg_electrode_locations{end+1} = [curr_x;curr_y]; % keep track of locations for aligning electrodes.
    wg_electrode_total_length{end+1} = device.wg_length;

    device_counter = device_counter + 1;

end
curr_x = curr_x + P.wg_horz_spacing*0.8 - 14;
clear device;
master_pattern_subchip.addelement(wg_electrode_grp);

%% Make straight waveguide without electrode for optical alignment on other side of the chip
N_wgs = 1;
grating_holes_grp = Group([0,0],[1;0],[0;1],{});

%DEVICE PARAMETERS
wg_width = P.wg_width;
% reflector parameters
scaling_holes = repmat([1 0.97 1.03 0.94 1.05], 1, 1);

device = {};

for ii = 1:N_wgs
    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.width_corr = true;
    device.wg_length = P.total_wg_length;
    device.tag_text = num2str(device_counter);
    device.scaling_holes = scaling_holes(ii);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,device.scaling_holes,P.wg_width,P.si_protect_width);
    
    wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    grating_holes_grp.addelement(wg_dev);

    curr_x = curr_x + P.write_field_x;
    device_counter = device_counter + 1;

end
%curr_x = curr_x + P.array_distance;
clear device;
master_pattern_subchip.addelement(grating_holes_grp);

%% Make fishbones w/out electrodes
N_fishbones = 10;
fishbone_grp = Group([0,0],[1;0],[0;1],{});

%DEVICE PARAMETERS
wg_width = P.fishbone_wg_width;
cavity_length = 450*ones(1,N_fishbones);
% reflector parameters
N_ramp_in = 10 * ones(1, N_fishbones);
N_mirror = repmat([12 9 6 3 0], 1, 2);
N_ramp_out = 10 * ones(1, N_fishbones);
ampl = 1.36*[0.100*ones(1, 5) 0.090*ones(1,5)];
period = (lambda0/(2*neff))*ones(1,N_fishbones);

device = {};
curr_x = curr_x + P.write_field_x;

for ii = 1:N_fishbones
    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.cavity_length = cavity_length(ii);
    % reflector parameters
    device.N_ramp_in = N_ramp_in(ii);
    device.N_mirror = N_mirror(ii);
    device.N_ramp_out = N_ramp_out(ii);
    device.ampl = ampl(ii);
    device.period = period(ii);
    device.ramp_type = 'lin';
    device.tag_text = num2str(device_counter);
    device.type = 'rect';
    device.duty_cycle = 0.70;
    device.width_corr = true;

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,device.wg_width,P.si_protect_width);
    
    fishbone_dev = makeFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B);
    fishbone_grp.addelement(fishbone_dev);

    curr_x = curr_x + P.write_field_x;
    device_counter = device_counter + 1;

end
clear device;
master_pattern_subchip.addelement(fishbone_grp);


curr_x = curr_x + 2*P.write_field_x*0.7;

%% Make ring resonator
ring_grp = Group([0,0],[1;0],[0;1],{});

curr_x = curr_x + P.write_field_x + 17;

%Make a standard Si strip ring resonator
device = {};

N_rings = 1;
for ii = 1:N_rings

    device.wg_length = P.total_wg_length;
    device.wg_width = P.wg_width;
    device.ring_cpl_spacing = 0.300;
    device.ring_radius = P.ring_gc_horz_spacing/2 - device.wg_width - device.ring_cpl_spacing ;

    device.tag_text1 = num2str(device_counter);
    device.tag_text2 = num2str(device_counter + 1);

    %GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes);
    GC_fn = @() gratingcoupler_airSOI220_27deg_twostepmeta_focusing(L.SI_LAYER,L.HOLE_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER,B.corrections,B.hide_GC_holes,P.scaling_pitch,P.scaling_holes,device.wg_width,P.si_protect_width);
    
    
    ring_grp.add( makeDoubleWaveguideRingDevice(curr_x, curr_y, GC_fn, GC_fn, device, P, L, B) );
    device_counter = device_counter + 2;
    curr_x = curr_x + 2*P.write_field_x;
end
clear device

ring_grp.translate([-P.ring_gc_horz_spacing/2;0]); % translate ring group back to center ring in write field.
master_pattern_subchip.addelement(ring_grp);


%% Make chip dicing marks

dice_mark_grp = Group([0,0],[1;0],[0;1],{});

diced_chip_width = 7000;
diced_chip_height = 3000;

% find center location:
n = length(wg_electrode_locations)/2;
center = mean([wg_electrode_locations{n} wg_electrode_locations{n+1}]');

centerX = center(1);
centerY = origin(2);

P.minX = centerX - diced_chip_width/2;
P.maxX = centerX + diced_chip_width/2;
P.minY = origin(2) - diced_chip_height/2;
P.maxY = origin(2) + diced_chip_height/2;

align_inner_width = 2;
align_outer_width = 10;
align_inner_length = 4;
align_arm_length = 150; % = total length of one arm (half of full cross length).

if B.include_dicing_marks

    % make PQRS marks (nested design).
    dice_mark_grp.addelement( NestedAlignmentMark([P.minX ; centerY], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
    dice_mark_grp.addelement( NestedAlignmentMark([P.maxX ; centerY], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
    dice_mark_grp.addelement( NestedAlignmentMark([centerX; P.maxY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
    dice_mark_grp.addelement( NestedAlignmentMark([centerX; P.minY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
    dice_mark_grp.addelement( NestedAlignmentMark([P.maxX ; P.minY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
    dice_mark_grp.addelement( NestedAlignmentMark([P.minX ; P.maxY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
    dice_mark_grp.addelement( NestedAlignmentMark([P.minX ; P.minY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
    dice_mark_grp.addelement( NestedAlignmentMark([P.maxX ; P.maxY ], L.DICE_MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

    master_pattern_subchip.addelement(dice_mark_grp);
end

%% Make large pads and electrodes for wirebonding

electrode_wires_grp = Group([0,0],[1;0],[0;1],{});
rail_offset = 800;
dicing_offset = 200;
P.inner_minX = P.minX + dicing_offset;
P.inner_maxX = P.maxX - dicing_offset;
P.inner_minY = P.minY + dicing_offset;
P.inner_maxY = P.maxY - dicing_offset;

P.big_pad_height = 400;
P.big_pad_width = 500;
P.rail_width = 100;
P.metal_wire_width = 10;
P.alignment_slop = 15;
%P.electrode_offset = electrode_gap;
pad_horz_offset = P.inner_maxX - 300 - centerX;
pad_horz_compens = 375;
pad_vert_compens = -100;
P.ground_width = 5250;

% make contact pads
BL = [centerX - pad_horz_offset; -rail_offset - P.big_pad_height];
TR = [centerX - pad_horz_offset + P.big_pad_width; -rail_offset];

BL_pad = CornerRect(BL,TR,L.HB_METAL_POS_LAYER); % BL pad
BL_pad.translate(pad_horz_compens,pad_vert_compens);
electrode_wires_grp.add(BL_pad);

BR_pad = CornerRect(BL,TR,L.HB_METAL_POS_LAYER);
BR_pad.translate(pad_horz_compens,pad_vert_compens);
BR_pad.mirror([centerX; centerY], [0;1])
electrode_wires_grp.add(BR_pad);

BL_ground = [centerX - P.ground_width/2; -rail_offset - P.big_pad_height];
TR_ground = [centerX + P.ground_width/2; -rail_offset];

T_pad = CornerRect(BL_ground,TR_ground, L.HB_METAL_POS_LAYER);
T_pad.mirror([centerX; centerY], [1;0])
T_pad.translate(0,-pad_vert_compens);
electrode_wires_grp.add(T_pad);

if B.split_rails
    % make rails
    rail_gap = 150;
    TR = [centerX - rail_gap/2; -rail_offset];
    BL = [centerX - pad_horz_offset + P.big_pad_width; -rail_offset - P.rail_width];

    BL_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER); % BL rail
    electrode_wires_grp.add( BL_rail);

    BR_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER);
    BR_rail.mirror([centerX; centerY], [0;1])
    electrode_wires_grp.add( BR_rail);

    TL_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER);
    TL_rail.mirror([centerX; centerY], [1;0])
    electrode_wires_grp.add( TL_rail);

    TR_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER);
    TR_rail.mirror([centerX; centerY], [1;0])
    TR_rail.mirror([centerX; centerY], [0;1])
    electrode_wires_grp.add( TR_rail);

else
    % make rails
    TR = [centerX + pad_horz_offset - P.big_pad_width - 375; -rail_offset];
    BL = [centerX - pad_horz_offset + P.big_pad_width + 375; -rail_offset - P.rail_width];

    bottom_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER);
    bottom_rail.translate(0,pad_vert_compens);
    electrode_wires_grp.add( bottom_rail);

%     top_rail = CornerRect(BL,TR, L.HB_METAL_POS_LAYER);
%     top_rail.mirror([centerX; centerY], [1;0])
%     electrode_wires_grp.add( top_rail);

end

%% Make small electrodes next to the waveguides

for ii = 1:length(wg_electrode_locations)

    curr = wg_electrode_locations{ii};
    curr_x = curr(1);
    curr_y = curr(2);
    
    % params
    wg_width = wg_widths(ii);

    % % make EOS6 electrode part
    EOS6_elec_wire_width = 10;
    EOS6_elec_wire_length = 150; % this is the horizontal electrical wire length

    % define some points (to simplify making electrode structure) (making
    % right side of electrode)
    Ax1 = wg_width/2 + P.electrode_offset;
    Ay1 = wg_electrode_total_length{ii}/2;
    Ax2 = Ax1 + EOS6_elec_wire_width;
    Ay2 = EOS6_elec_wire_width/2;
    Ax3 = EOS6_elec_wire_length - EOS6_elec_wire_width/2;
    Ax4 = EOS6_elec_wire_length + EOS6_elec_wire_width/2;

    points = {[Ax1;Ay1], [Ax2;Ay1], [Ax2;Ay2], [Ax3;Ay2], [Ax3;Ay1], [Ax4;Ay1], [Ax4;-Ay2], [Ax2;-Ay2], [Ax2;-Ay1], [Ax1;-Ay1]};

    TR_EOS6_electrode = Polygon(points, L.EOS6_METAL_LAYER);
    TR_EOS6_electrode.translate([curr_x;curr_y]);
    electrode_wires_grp.add(TR_EOS6_electrode);
    BR_EOS6_electrode = Polygon(points, L.EOS6_METAL_LAYER);
    BR_EOS6_electrode.translate([curr_x;curr_y]);
    BR_EOS6_electrode.rotate([curr_x;curr_y], 180);
    electrode_wires_grp.add(BR_EOS6_electrode);


    % % make EOS3 electrode part
    EOS3_elec_wire_width = 40;

    TR = [Ax3 + EOS3_elec_wire_width/2; rail_offset];
    BL = [Ax3 - EOS3_elec_wire_width/2; -Ay2];


    TR_EOS3_electrode = CornerRect(BL,[TR(1);TR(2) - pad_vert_compens], L.HB_METAL_POS_LAYER);
    TR_EOS3_electrode.translate([curr_x;curr_y]);
    electrode_wires_grp.add(TR_EOS3_electrode);
    BL_EOS3_electrode = CornerRect(BL,[TR(1);TR(2) - pad_vert_compens], L.HB_METAL_POS_LAYER);
    BL_EOS3_electrode.translate([curr_x;curr_y]);
    BL_EOS3_electrode.rotate([curr_x;curr_y], 180);
    electrode_wires_grp.add(BL_EOS3_electrode);

    %     % define some points
    %     Bx1 = Ax3 - EOS3_elec_wire_width;
    %     Bx2 = Ax3;
    %     Bx3 = Bx1 + EOS3_pad_width;
    %     By1 = -Ay2;
    %     By3 = EOS3_elec_wire_length;
    %     By2 = By3 - EOS3_pad_width;
    %
    %     points = {[Bx1;By1], [Bx1;By3], [Bx3;By3], [Bx3;By2], [Bx2;By2], [Bx2;By1]};
    %
    %     TR_EOS3_electrode = Polygon(points, L.HB_METAL_POS_LAYER);
    %     TR_EOS3_electrode.translate([curr_x;curr_y]);
    %     electrode_wires_grp.add(TR_EOS3_electrode);
    %     BR_EOS3_electrode = Polygon(points, L.HB_METAL_POS_LAYER);
    %     BR_EOS3_electrode.translate([curr_x;curr_y]);
    %     BR_EOS3_electrode.rotate([curr_x;curr_y], 180);
    %     electrode_wires_grp.add(BR_EOS3_electrode);


end

master_pattern_subchip.addelement(electrode_wires_grp);


%% Make ebeam/heidelberg alignment marks

if B.include_align_marks
    centerX = (P.inner_maxX + P.inner_minX) / 2;
    centerY = origin(2);
    align_mark_grp = makeProtectedAlignmentMarks_v2(centerX, centerY, P, L);
    master_pattern_subchip.add(align_mark_grp);
end

%% Make Stanford logo

if B.include_logos
    logo_grp = Group([0,0],[1;0],[0;1],{}, L.HB_SILICON_PROTECT_LAYER);

    boxX = 1.5;  % box size of 1 um gives logos which are roughly 150 um x 200 um.  Origin is in top left of logo.
    boxY = 1.5;
    X = 0;
    Y = 0;

    my_tree = Group([0;0], [1;0], [0;1], process_tree(boxX,boxY,X,Y));
    my_tree.translate([P.inner_minX + 1000 + 2870; P.maxY - 200 - 2348]);
    logo_grp.add(my_tree);

    logo_scl_factor = 0.86; % make logos the same size
    my_logo = Group([0;0], [1;0], [0;1], process_logo(boxX*logo_scl_factor,boxY*logo_scl_factor,X,Y));
    my_logo.translate([P.inner_minX + 1600 + 2870; P.maxY - 200 - 2348]);
    logo_grp.add(my_logo);

    master_pattern_subchip.add(logo_grp);

    % make chip tag
    tag_pos = [P.inner_maxX - 1000 + 520; P.maxY - 300] ;
    logo_grp.add(makeTag(tag_pos, chip_tag, B, L.HB_SILICON_PROTECT_LAYER, 'box_size', 10));
end

%% Make rectangle bounding box for metal protection

master_pattern_subchip.add( Rect([centerX;centerY], 7000, 2900, L.HB_METAL_GAP_LAYER) );

%% Re-center chip at (0,0)

master_pattern_subchip.translate(-centerX, -centerY);

%% WRITE DFX FILE
if B.generate_CAD
    display('Writing layers')
    master_pattern_subchip.todxf(filename);
    % files not closing properly, so close all here.
    fclose('all')
end

out = master_pattern_subchip;

toc

end


%%%%%%%%%%%%%% Functions start here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function align_mark_grp = makeProtectedAlignmentMarks(centerX, centerY, P, L)

align_mark_grp = Group([0,0],[1;0],[0;1],{});

align_inner_width = 0.1;
align_outer_width = 5;
align_inner_length = 4;
align_arm_length = 100; % = total length of one arm (half of full cross length).

% make PQRS + corner marks (nested design).
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

% add two layer align mark structures to determine alignment between
% EOS6 metal and EOS6 photonics write.
align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );
align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX + 2*align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );


% make a copy of all these alignment marks, but protect them from Si etch
align_inner_width = 20;
align_outer_width = 20;
align_inner_length = 4;
align_arm_length = 100;

align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

align_mark_grp.addelement( Rect([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length],  align_inner_width, align_inner_width, L.HB_SILICON_PROTECT_LAYER) );
align_mark_grp.addelement( Rect([centerX + 2*align_arm_length; P.inner_minY + align_arm_length],  align_inner_width, align_inner_width, L.HB_SILICON_PROTECT_LAYER) );


end

function align_mark_grp = makeProtectedAlignmentMarks_v2(centerX, centerY, P, L)

align_mark_grp = Group([0,0],[1;0],[0;1],{});

align_inner_width = 0.1;
align_outer_width = 5;
align_inner_length = 4;
align_arm_length = 100; % = total length of one arm (half of full cross length).

% make PQRS + corner marks (nested design).
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
%align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

% add two layer align mark structures to determine alignment between
% EOS6 metal and EOS6 photonics write.
%align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length], L.NEG_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );
align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX + 2*align_arm_length; P.inner_minY + align_arm_length], L.NEG_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );


% make a copy of all these alignment marks, but protect them from Si etch
align_inner_width = 20;
align_outer_width = 20;
align_inner_length = 4;
align_arm_length = 100;

align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
%align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
%align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
%align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.HB_UNDERCUT, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR


%align_mark_grp.addelement( Rect([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length],  align_inner_width, align_inner_width, L.HB_SILICON_PROTECT_LAYER) );
align_mark_grp.addelement( Rect([centerX + 2*align_arm_length; P.inner_minY + align_arm_length],  align_inner_width, align_inner_width, L.HB_SILICON_PROTECT_LAYER) );


end


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
      elements{end+1} = makeTag_v2(tag_pos, device.tag_text, B, L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
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
      elements{end+1} = makeTag_v2(tag_pos, device.tag_text, B, L.NEG_LAYER);
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
      elements{end+1} = makeTag_v2(tag_pos, device.tag_text, B, L.NEG_LAYER);
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
      elements{end+1} = Rect([curr_x + device.wg_width + device.ring_cpl_spacing + device.ring_radius; curr_y], 2*r, 2*r, L.HB_SILICON_PROTECT_LAYER);
  end

  %make left waveguide
  device.tag_text = device.tag_text1;
  wg_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn1, device, P, L, B);
  elements{end+1} = wg_dev;

  %make ring resonator
  elements{end+1} = RingRes([curr_x + device.wg_width + device.ring_cpl_spacing + device.ring_radius; curr_y], device.ring_radius, device.wg_width, P.trench, P.trench, P.trench, P.N_points, L.SI_LAYER, L.NEG_LAYER, L.HOLE_LAYER, L.INNER_POS_LAYER);

  if B.include_label
      tag_pos = [curr_x; curr_y  - device.wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag_v2(tag_pos, device.tag_text1, B, L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
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
      elements{end+1} = Rect([0; 0], width, device.wg_length, L.HB_SILICON_PROTECT_LAYER);
  end

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - device.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - device.wg_length/2], -90);

  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x; curr_y - device.wg_length/2]);
  elements{end}.rotate([curr_x; curr_y - device.wg_length/2], -90);
  %elements{end}.mirror([curr_x;curr_y], [1;0]); % MIRROR function shifts
  %the holes in weird way!!!
  elements{end}.rotate([curr_x;curr_y],180);

  elements{end+1} = makeStripWaveguide(device.wg_width, device.wg_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  if B.include_label
      tag_pos = [curr_x; curr_y  - device.wg_length/2 - P.taper_length ] + P.tag_offset ;
      elements{end+1} = makeTag_v2(tag_pos, device.tag_text, B, L.NEG_LAYER, L.HB_SILICON_PROTECT_LAYER);
      %elements{end+1} = Rect(tag_pos, 30, 15, L.HB_SILICON_PROTECT_LAYER);
  end
  % make a Group out of the elements
  thru_wg = Group([x0,y0],[1;0],[0;1],elements);
end
