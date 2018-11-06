% QuarterWaveQEOM_v1 Layout
% JDW 2017-09-26
% EOS6 photonics + Heidelberg + Heidelberg + EOS6 metal

function out = QuarterWaveQEOM_v1_subchip(B, L_coupler, chip_tag)

addpath(genpath('C:\Users\jwitmer\Documents\GitHub\CAD'))
addpath(genpath('Z:\jwitmer\CAD_working_folder\2017_09_26_QuarterWaveQEOM_v1\McKenna_functions'))
% addpath(genpath('/home/tim/github_local/CAD'))
% addpath(genpath('/media/user_data/jwitmer/CAD_working_folder/2017_09_26_QuarterWaveQEOM_v1/'))

%cd('Z:\jwitmer\CAD_working_folder\2017_09_26_QuarterWaveQEOM_v1')

% clear
tic

% dxf file name, without .dxf
filename = 'QuarterWaveQEOM_v1';

%% LAYERS
% L for Layers
L.SI_LAYER = 'SiLayer';
L.NEG_LAYER = 'NegLayer';
L.HOLE_LAYER = 'HoleLayer'; % Si to be removed from positive regions.
L.INNER_POS_LAYER = 'InnerPosLayer'; % for ring inside

L.HB_SI_PROTECT_LAYER = 'HBSiProtectLayer';
L.HB_METAL_GAP_LAYER = 'HBMetalGapLayer';
L.HB_METAL_POS_LAYER = 'HBMetalPosLayer';
L.EOS6_METAL_LAYER = 'EOS6MetalLayer';

L.LOGO_LAYER = L.HB_METAL_GAP_LAYER;

%% initialize master group
master_group    =  Group([0,0],[1;0],[0;1],{});

%% Parameters:

elec_spacing = 2.7;
elec_wire_width = 3;
elec_height = 450;
meander_X_offset = 400;



%% Define converter center positions

converter_x_pos = linspace(-2580, 1740, 5);
converter_y_pos = 0;

res_length = linspace(5000, 5800, 5);


for ii = 1:length(converter_x_pos)

    %% Make converter cut-out box
    cutout_h = 1150;
    cutout_w = 300;
    cutout_vert_shift = 75;
    master_group.add( Rect( [converter_x_pos(ii); converter_y_pos + cutout_vert_shift], cutout_w, cutout_h, L.HB_METAL_GAP_LAYER));

    %% make iith optical device here, with cavity centered at [converter_x_pos(ii); converter_y_pos]   %%%%%%%%%%%%

    if chip_tag ~= 'C'
      if ii ~= 5
        photon_dev = makePhotonDevice(ii, chip_tag, B.hide_GC_holes);
        photon_dev.translate(converter_x_pos(ii), converter_y_pos);
        master_group.add(photon_dev);
      end
    end
    %% Make meander resonators
    % c : top left corner of bounding box
    % s : span (width) of bounding box
    % l : total length of meander
    % w : width of wire
    % l_c : length of coupler segment
    % r : turn radius for elbow segments
    % N : number of points to defined arc paths
    % gtw : gap-to-wire ratio

    c = [-1000; 1000];
    res_L = res_length(ii) - (meander_X_offset - cutout_w/2) ;
    w = 6.4; % should come out to 5 um;
    gap_w = 11.7; % should come out to 13 um
    l_c = L_coupler;
    r = 30;
    N = 50;
    gtw = gap_w/w;
    N_turns = 12;

    % set s so that meander ends precisely at center of winding:
    % l = l_c + N_turns*( (s-2*r) + pi*r) + pi*r + s/2;

    s = (res_L - l_c + (2*N_turns - N_turns*pi - pi + 1)*r) / (N_turns + 1/2);
    % s = 10000;

    next_meander = SimpleMeanderGaps_01(c, s, res_L, w, l_c, r, N, gtw, L.HB_METAL_GAP_LAYER);
    next_meander.mirror([0;0],[1;0]);
    ext = next_meander.getExtent;
    BL = ext{1};
    TL = ext{4};
    meander_height = TL(2) - BL(2);

    % bring meander to correct position
    next_meander.translate( -TL + [converter_x_pos(ii) + meander_X_offset - s/2 - w/2 - w*gtw; w/2 + gtw*w]  );

    master_group.add( next_meander );

    % add gaps extension to converter box:
    gaps_path = {[converter_x_pos(ii) + cutout_w/2; converter_y_pos], [converter_x_pos(ii) + meander_X_offset; converter_y_pos]};
    master_group.add( Gaps([0,0], gaps_path, w, gtw, L.HB_METAL_GAP_LAYER) );

    % make right and left electrodes:
    electrode_input_length = meander_X_offset*0.8;
    right_electrode = TElectrode([converter_x_pos(ii) + elec_spacing/2; converter_y_pos], elec_height, electrode_input_length, elec_wire_width, w, L.EOS6_METAL_LAYER);
    left_electrode = TElectrode([converter_x_pos(ii) + elec_spacing/2; converter_y_pos], elec_height, electrode_input_length, elec_wire_width, w, L.EOS6_METAL_LAYER);
    left_electrode.mirror([converter_x_pos(ii); converter_y_pos], [0;1]);

    master_group.add(right_electrode);
    master_group.add(left_electrode);

    % make large wirebond capacitor
    bond_pad_height = cutout_h;
    bond_pad_width = 125;
    bond_pad_gap = 5;
    bond_pad_input_wire_length = 80;
    bond_pad_input_wire_width = w;
    bond_pad_vert_offset = cutout_vert_shift;
    master_group.add( SquareBondPadCap([-cutout_w/2 + converter_x_pos(ii); converter_y_pos],...
        bond_pad_height, bond_pad_width, bond_pad_gap, bond_pad_input_wire_length,...
        bond_pad_input_wire_width, bond_pad_vert_offset,  L.HB_METAL_GAP_LAYER, L.HB_METAL_POS_LAYER) );


end


%% Make feedline

% l_1 : length of wide launcher segment
% l_2 : taper length
% L1 : length of vertical segment
% L2 : length of horizontal segment
% R : bend radius
% N : number of points for each arc segment
% W : width of wide segment of launcher
% w : width of feedline wire
% gtw : gap-to-wire ratio

l_1 = 200;
l_2 = 200;
L1 = 50;
R = 100;
N = 50;
W = 200;
L2 = 4920 - 2*R;
feed_w = 21.4;  % should come out to 20 um
feed_gap_w = 4.7; % should come out to 6 um
feed_gtw = feed_gap_w / feed_w;

feedline_y = converter_y_pos + w/2 + w*gtw - meander_height - feed_w/2 - feed_w*feed_gtw + w*gtw;

feed_start = [-L2/2 - R; feedline_y - l_1 - l_2 - L1 - R];

feedline = UFeedlineTwoLayer(feed_start, l_1, l_2, L1, L2, R, N, W, feed_w, feed_gtw, N, L.HB_METAL_GAP_LAYER, L.HB_METAL_POS_LAYER);
master_group.add( feedline );

feedline_ext = feedline.getExtent;
feedline_BL = feedline_ext{1};


%% Make bounding box
box_corner = [-3500; feedline_BL(2) - 150];
box_l = 6700;
box_w = 2600;
box_L = 7000;
box_W = 2900;

bounding_box = BoundingBox_01(box_corner, box_w, box_l, box_W, box_L);
bounding_box.layer = L.HB_METAL_GAP_LAYER;
master_group.add( bounding_box );




%%  if desired, make extra optical devices here (ring resonators, thru wg's, etc) %%%%%%%%%%%%

if chip_tag ~= 'C'
  photon_dev_group = makeSidePhotonDevices(chip_tag, B.hide_GC_holes);
  photon_dev_group.translate(2580, 0);
  master_group.add(photon_dev_group);
end


%% Add alignment marks

if B.include_align_marks

    ext = bounding_box.getExtent;
    BL = ext{1};
    TR = ext{3};

    centerX = (TR(1) + BL(1) )/2;
    centerY = (TR(2) + BL(2) )/2;

    align_mark_offset = 200; % offset of alignment mark edge from chip border edge

    P.inner_minX = BL(1) + align_mark_offset;
    P.inner_minY = BL(2) + align_mark_offset;
    P.inner_maxX = TR(1) - align_mark_offset;
    P.inner_maxY = TR(2) - align_mark_offset;

    align_mark_grp = makeProtectedAlignmentMarks(centerX, centerY, P, L);

    master_group.add(align_mark_grp);

end

%%

if B.include_logos

    logo_grp = Group([0,0],[1;0],[0;1],{}, L.LOGO_LAYER);

    boxX = 1.5;  % box size of 1 um gives logos which are roughly 150 um x 200 um.  Origin is in top left of logo.
    boxY = 1.5;
    X = 0;
    Y = 0;

    % make Stanford logo
    my_tree = Group([0;0], [1;0], [0;1], process_tree(boxX,boxY,X,Y));
    my_tree.translate([P.inner_minX + 350; P.inner_maxY - 40]);
    logo_grp.add(my_tree);

    % make LINQS logo
%     logo_scl_factor = 0.86; % make logos the same size
%     my_logo = Group([0;0], [1;0], [0;1], process_logo(boxX*logo_scl_factor,boxY*logo_scl_factor,X,Y));
%     my_logo.translate([P.inner_minX + 1600; P.inner_maxY - 200]);
%     logo_grp.add(my_logo);

    % make chip tag
    tag_pos = [P.inner_maxX - 500; P.inner_maxY - 60] ;
    logo_grp.add(makeTag(tag_pos, chip_tag, B, L.LOGO_LAYER, 'box_size', 10));

    master_group.add(logo_grp);

end

%% Re-center chip at (0,0)

% master_group.translate(-centerX, -centerY);


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

align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; centerY], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; centerY], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_minY + align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_maxY - align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length; P.inner_minY + align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length; P.inner_maxY - align_arm_length], L.HB_SI_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

align_mark_grp.addelement( Rect([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length],  align_inner_width, align_inner_width, L.HB_SI_PROTECT_LAYER) );
align_mark_grp.addelement( Rect([centerX + 2*align_arm_length; P.inner_minY + align_arm_length],  align_inner_width, align_inner_width, L.HB_SI_PROTECT_LAYER) );


end
