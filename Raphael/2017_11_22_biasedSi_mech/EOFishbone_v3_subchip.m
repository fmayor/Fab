% EOFishbone_v3 Layout
% JDW 2017-08-10
% EOS6 + EOS 3

function out = EOFishbone_v3_subchip(B, electrode_gap, chip_tag)

addpath(genpath('C:\Users\jwitmer\Documents\GitHub\CAD'))
addpath(genpath('Z:\jwitmer\CAD_working_folder\2017_09_26_QuarterWaveQEOM_v1\McKenna_functions'))
% cd('Z:\jwitmer\CAD_working_folder\2017_08_10_EO_fishbone_v3')

% clear
tic

% dxf file name, without .dxf
filename = '/Users/Felix/Documents/Master Thesis/Fab/Electromechanical_fin/CAD/generation1/EOFishbone_v3_release';

B.heidelberg_si_protect = true;
P.si_protect_width = 1.5;

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
%L.EOS6_METAL_LAYER = 'HBMetalPosLayer';

L.MODULATION_ELEC_EOS3_LAYER = L.HB_METAL_POS_LAYER;
L.MODULATION_ELEC_EOS6_LAYER = L.EOS6_METAL_LAYER;
L.MISC_METAL_LAYER = L.EOS6_METAL_LAYER;
L.DICE_MARK_LAYER = 'DiceMarkLayer';

L.LOGO_LAYER = L.HB_SILICON_PROTECT_LAYER;




%% GLOBAL PARAMETERS
% P for Parameters
P.write_field_x = 60;
P.write_field_y = 60;


P.bragg_point_res = 0.010; % use fine resolution for CAD only, not for plotting


P.trench = 3; % trench around normal waveguides
P.taper_length = 345; % for grating couplers
P.total_wg_length = 270 + 290; % length between end of GC tapers (default: 400)
P.scaling_pitch = 1.0;  %dev16 from v1 chip
P.tag_offset = [-10;-40];
P.wg_horz_spacing = 480; %P.write_field_x;  %space waveguides by the write field
P.fishbone_wg_width = 0.372;
P.ring_gc_horz_spacing = 46;
P.wg_width = 0.476; % corrected, should come out around 450 nm.
P.N_points = 3;  %this is for the RingRes function. NEED TO REMOVE FROM FNC!
P.contact_wire_width = 1; % width of wire protruding from slot mode converter
P.contact_wire_length = 80;
P.GC_hole_corr = false;

P.fin_width = 0.100*1.07;
P.fin_gap_fishbone=0.100; %distance fin to fishbone
P.fin_gap_electrode=[0.750,0.250,-0.250,1.0]; %distance fin to left electrode
P.fin_length=450-2*11; %device.cavity_length - space for tapers
P.taper_length = 9; %quarter period of cosine
P.taper_amplitude = 41.6; %needs to be adjusted if change fin_width or taper_length

%% Initialize
master_group    =  Group([0,0],[1;0],[0;1],{});

% pattern must be away from [0,0] to avoid errors in AutoCAD, so offset
% origin
origin = [0;0];

curr_x = origin(1) + P.wg_horz_spacing/2;
curr_y = origin(2);

device_counter = 1;
fishbone_locations = {};
fishbone_total_length = []; % length of cavity, including reflectors, for each fishbone device.


%% make the fishbone array

N_fishbones = 4;
fishbone_grp = Group([0,0],[1;0],[0;1],{});

% want to shift center wavelength from 1540 to 1578, so new period is 434
% nm


%DEVICE PARAMETERS
wg_width = P.fishbone_wg_width;
cavity_length = 450*ones(1, N_fishbones);
% reflector parameters
N_ramp_in = 10 * ones(1, N_fishbones);
%N_mirror = [6 6 6 5 4 3 3 1]; %%%%%%
N_mirror = [3 3 3 3];
N_ramp_out = 10 * ones(1, N_fishbones);
ampl = 0.100 * ones(1, N_fishbones);
period = 0.436;

device = {};
curr_x = curr_x - P.wg_horz_spacing/2;
for ii = 1:N_fishbones

    curr_x = curr_x + P.wg_horz_spacing;
    if ii==2
       %curr_x=curr_x-electrode_gap-P.fishbone_wg_width/2+P.fin_gap_electrode(2)+P.fin_width+P.fin_gap_fishbone;
       curr_x=curr_x-electrode_gap-P.fishbone_wg_width/2;
    end

    %DEVICE PARAMETERS
    device.wg_width = wg_width;
    device.cavity_length = cavity_length(ii);
    % reflector parameters
    device.N_ramp_in = N_ramp_in(ii);
    device.N_mirror = N_mirror(ii);
    device.N_ramp_out = N_ramp_out(ii);
    device.ampl = ampl(ii);
    device.period = period;
    device.ramp_type = 'lin';
    device.tag_text = num2str(device_counter);
    device.type = 'rect';
    device.duty_cycle = 0.70;
    device.width_corr = true;

    GC_fn = @() makeTwoStepGrating([0; 0], device.wg_width, P.trench, P.taper_length, P.scaling_pitch, L.SI_LAYER, L.HOLE_LAYER, L.NEG_LAYER, P.GC_hole_corr, B.hide_GC_holes, B.heidelberg_si_protect, L.HB_SILICON_PROTECT_LAYER);

    fishbone_dev = makeFishboneCavityDevice(curr_x, curr_y, GC_fn, device, P, L, B);
%     if ii == 2
%       %make thru waveguide for RRRrrraphael
%       device.wg_width = 0.560;  %should around 500 nm
%       device.wg_length = P.total_wg_length;
%       fishbone_dev = makeStripWaveguideDevice(curr_x, curr_y, GC_fn, device, P, L, B);
%     end
    fishbone_grp.addelement(fishbone_dev);

    fishbone_locations{end+1} = [curr_x;curr_y]; % keep track of locations for aligning electrodes.
    fishbone_total_length(end+1) = device.cavity_length + 2*device.period*(device.N_ramp_in + device.N_mirror + device.N_ramp_out);

    device_counter = device_counter + 1;
    if ii==2
       %curr_x=curr_x+electrode_gap+P.fishbone_wg_width/2-P.fin_gap_electrode(2)-P.fin_width-P.fin_gap_fishbone;
       curr_x=curr_x+electrode_gap+P.fishbone_wg_width/2;
    end
end
curr_x = curr_x + P.wg_horz_spacing/2;
clear device;
master_group.addelement(fishbone_grp);

%% Make the fins
%need rectangle and smooth tapers
%curr_x = curr_x - P.wg_horz_spacing/2; 
fishbone_locations2=cell2mat(fishbone_locations);

fin_grp = Group([0,0],[1;0],[0;1],{});
fin_grp.layer=L.SI_LAYER;

for ii=1:N_fishbones
    %Fin:
    curr_x = fishbone_locations2(1,ii);
    curr_y = fishbone_locations2(2,ii);
    fin=Rect(curr_x-P.fin_gap_fishbone-P.fin_width-P.fishbone_wg_width/2,curr_y-P.fin_length/2,P.fin_width,P.fin_length); %position x, position y, width, height
    fin_grp.addelement(fin);

    %Tapers:
    curr_x=curr_x-P.fin_gap_fishbone-P.fin_width-P.fishbone_wg_width/2;
    curr_y=curr_y+P.fin_length/2;
    y_i=(curr_y:0.004:curr_y+P.taper_length);
    x_i=curr_x-P.taper_amplitude*(1-cos((y_i-curr_y)./(4*P.taper_length))); %taper function is a quarter period of a cosine of period C

    points=[x_i(1);y_i(1)];
    for ii=2:length(y_i)
        points=horzcat(points,[x_i(ii);y_i(ii)]);
    end

    for ii=1:length(y_i)
        points=horzcat(points,[x_i(length(y_i)+1-ii)+P.fin_width;y_i(length(y_i)+1-ii)]);
    end
    points=horzcat(points,[x_i(1);y_i(1)]);
    taperup=Polygon(points); %upper taper

    points2=vertcat(points(1,:),-points(2,:));
    taperdown=Polygon(points2); %lower taper

    fin_grp.addelement(taperup);
    fin_grp.addelement(taperdown);

end


master_group.addelement(fin_grp);
%% Make chip dicing marks


dice_mark_grp = Group([0,0],[1;0],[0;1],{});

diced_chip_width = 7000;
diced_chip_height = 3000;

% find center location:
n = length(fishbone_locations)/2;
center = mean( [fishbone_locations{n} fishbone_locations{n+1}]');

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

    master_group.addelement(dice_mark_grp);


end


%% Make electrical contacts for poling.

electrode_wires_grp = Group([0,0],[1;0],[0;1],{});

dicing_offset = 200;

P.inner_minX = P.minX + dicing_offset;
P.inner_maxX = P.maxX - dicing_offset;
P.inner_minY = P.minY + dicing_offset;
P.inner_maxY = P.maxY - dicing_offset;

rail_offset = 800;

P.big_pad_height = 300;
P.big_pad_width = 300;

P.rail_width = 50;
P.metal_wire_width = 10;
P.alignment_slop = 15;

P.electrode_offset = electrode_gap;

pad_horz_offset = P.inner_maxX - 300 - centerX;

% make contact pads
BL = [centerX - pad_horz_offset; -rail_offset - P.big_pad_height];
TR = [centerX - pad_horz_offset + P.big_pad_width; -rail_offset];

BL_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER); % BL pad
electrode_wires_grp.add( BL_pad);

BR_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
BR_pad.mirror([centerX; centerY], [0;1])
electrode_wires_grp.add( BR_pad);

TL_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
TL_pad.mirror([centerX; centerY], [1;0])
electrode_wires_grp.add( TL_pad);

TR_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
TR_pad.mirror([centerX; centerY], [0;1])
TR_pad.mirror([centerX; centerY], [1;0])
electrode_wires_grp.add( TR_pad);


if B.split_rails
    % make rails
    rail_gap = 150;
    TR = [centerX - rail_gap/2; -rail_offset];
    BL = [centerX - pad_horz_offset + P.big_pad_width; -rail_offset - P.rail_width];

    BL_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER); % BL rail
    electrode_wires_grp.add( BL_rail);

    BR_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    BR_rail.mirror([centerX; centerY], [0;1])
    electrode_wires_grp.add( BR_rail);

    TL_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    TL_rail.mirror([centerX; centerY], [1;0])
    electrode_wires_grp.add( TL_rail);

    TR_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    TR_rail.mirror([centerX; centerY], [1;0])
    TR_rail.mirror([centerX; centerY], [0;1])
    electrode_wires_grp.add( TR_rail);

else
    % make rails
    TR = [centerX + pad_horz_offset - P.big_pad_width; -rail_offset];
    BL = [centerX - pad_horz_offset + P.big_pad_width; -rail_offset - P.rail_width];

    bottom_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    electrode_wires_grp.add( bottom_rail);

    top_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    top_rail.mirror([centerX; centerY], [1;0])
    electrode_wires_grp.add( top_rail);

end

% make modulation electrodes.

for ii = 1:length(fishbone_locations)

    curr = fishbone_locations{ii};
    curr_x = curr(1);
    curr_y = curr(2);

    % % make EOS6 electrode part
    EOS6_elec_wire_width = 3; %5;
    EOS6_elec_wire_length = 150;

    % define some points (to simplify making electrode structure) (making
    % right side of electrode)
    Ax1 = P.fishbone_wg_width/2 + P.electrode_offset;
    % MODIFICATION OF ELECTRODES FOR RAPHAEL
    % DEV 7 AND 8 HAVE 3 um offset
    if ii == 7
      Ax1 = P.fishbone_wg_width/2 + 3.0;
    elseif ii == 8
      Ax1 = Ax1 + 0.394 + 1.5;
    end
    Ay1 = fishbone_total_length(ii)/2;
    Ax2 = Ax1 + EOS6_elec_wire_width;
    Ay2 = EOS6_elec_wire_width/2;
    Ax3 = EOS6_elec_wire_length - EOS6_elec_wire_width/2;
    Ax4 = EOS6_elec_wire_length + EOS6_elec_wire_width/2;

    points = {[Ax1;Ay1], [Ax2;Ay1], [Ax2;Ay2], [Ax3;Ay2], [Ax3;Ay1], [Ax4;Ay1], [Ax4;-Ay2], [Ax2;-Ay2], [Ax2;-Ay1], [Ax1;-Ay1]};

    TR_EOS6_electrode = Polygon(points, L.MODULATION_ELEC_EOS6_LAYER);
    TR_EOS6_electrode.translate([curr_x;curr_y]);
    electrode_wires_grp.add(TR_EOS6_electrode);
    BR_EOS6_electrode = Polygon(points, L.MODULATION_ELEC_EOS6_LAYER);
    BR_EOS6_electrode.translate([curr_x-electrode_gap+P.fin_gap_electrode(ii)+P.fin_width+P.fin_gap_fishbone;curr_y]);
    BR_EOS6_electrode.rotate([curr_x;curr_y], 180);
    electrode_wires_grp.add(BR_EOS6_electrode);


    % % make EOS3 electrode part
    EOS3_elec_wire_width = 20;

    TR = [Ax3 + EOS3_elec_wire_width/2; rail_offset];
    BL = [Ax3 - EOS3_elec_wire_width/2; -Ay2];


    TR_EOS3_electrode = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    TR_EOS3_electrode.translate([curr_x;curr_y]);
    electrode_wires_grp.add(TR_EOS3_electrode);
    BL_EOS3_electrode = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
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
    %     TR_EOS3_electrode = Polygon(points, L.MODULATION_ELEC_EOS3_LAYER);
    %     TR_EOS3_electrode.translate([curr_x;curr_y]);
    %     electrode_wires_grp.add(TR_EOS3_electrode);
    %     BR_EOS3_electrode = Polygon(points, L.MODULATION_ELEC_EOS3_LAYER);
    %     BR_EOS3_electrode.translate([curr_x;curr_y]);
    %     BR_EOS3_electrode.rotate([curr_x;curr_y], 180);
    %     electrode_wires_grp.add(BR_EOS3_electrode);


end


master_group.addelement(electrode_wires_grp);


%% Add alignment marks
% 
% if B.include_align_marks
% 
%     align_mark_grp = Group([0,0],[1;0],[0;1],{});
% 
%     align_inner_width = 0.1;
%     align_outer_width = 5;
%     align_inner_length = 4;
%     align_arm_length = 100; % = total length of one arm (half of full cross length).
% 
%     centerX = (P.inner_maxX + P.inner_minX) / 2;
%     centerY = origin(2);
% 
%     align_mark_grp = makeProtectedAlignmentMarks(centerX, centerY, P, L);
% 
%     master_group.add(align_mark_grp);
% 
% end

%%

% if B.include_logos
% 
%     logo_grp = Group([0,0],[1;0],[0;1],{}, L.LOGO_LAYER);
% 
%     boxX = 1.5;  % box size of 1 um gives logos which are roughly 150 um x 200 um.  Origin is in top left of logo.
%     boxY = 1.5;
%     X = 0;
%     Y = 0;
% 
%     my_tree = Group([0;0], [1;0], [0;1], process_tree(boxX,boxY,X,Y));
%     my_tree.translate([P.inner_minX + 1000; P.maxY - 200]);
%     logo_grp.add(my_tree);
% 
%     logo_scl_factor = 0.86; % make logos the same size
%     my_logo = Group([0;0], [1;0], [0;1], process_logo(boxX*logo_scl_factor,boxY*logo_scl_factor,X,Y));
%     my_logo.translate([P.inner_minX + 1600; P.maxY - 200]);
%     logo_grp.add(my_logo);
% 
%     master_group.add(logo_grp);
% 
% 
%     % make chip tag
%     tag_pos = [P.inner_maxX - 1000; P.maxY - 300] ;
%     logo_grp.add(makeTag(tag_pos, chip_tag, B, L.LOGO_LAYER, 'box_size', 10));
% 
% 
% end

%% Make rectangle bounding box for metal protection

master_group.add( Rect([centerX;centerY], 7000, 2900, L.HB_METAL_GAP_LAYER) );



%% Re-center chip at (0,0)

master_group.translate(-centerX, -centerY);



%% WRITE DFX FILE
if B.generate_CAD
    display('Writing layers')
    master_group.todxf(filename);
    % files not closing properly, so close all here.
    fclose('all')
end

out = master_group;

toc

end


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
