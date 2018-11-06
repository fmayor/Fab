% Fin driving Layout with fin cavities
% FMM from JDW+RVL+CJS 2018-03-09
% EOS6 + EOS3

function out = Olympus(B, chip_tag)

% clear
tic

% dxf file name, without .dxf
filename = 'Olympus';

B.heidelberg_si_protect = true;
P.si_protect_width = 1.5;

%% LAYERS
% L for Layers
L.SI_LAYER = 'SiLayer';
L.NEG_LAYER = 'NegLayer';
L.HOLE_LAYER = 'HoleLayer'; % Si to be removed from positive regions.
L.INNER_POS_LAYER = 'InnerPosLayer'; % for ring inside
L.HB_SILICON_PROTECT_LAYER = 'HBSiProtectLayer';% L.METAL_ON_SI_LAYER = 'MetalOnSiLayer'; % Si to be removed from positive regions.
%L.BANDAGE_LAYER = 'BandageLayer';
L.MARK_LAYER='MarkLayer';
L.HB_METAL_GAP_LAYER = 'HBMetalGapLayer';
L.HB_METAL_POS_LAYER = 'HBMetalPosLayer';
L.EOS6_METAL_LAYER = 'SpiralLayer';
%L.EOS6_METAL_LAYER = 'HBMetalPosLayer';

L.MODULATION_ELEC_EOS3_LAYER = L.HB_METAL_POS_LAYER;
L.MODULATION_ELEC_EOS6_LAYER = L.EOS6_METAL_LAYER;
L.MISC_METAL_LAYER = L.EOS6_METAL_LAYER;
L.DICE_MARK_LAYER = 'DiceMarkLayer';

L.LOGO_LAYER = L.HB_SILICON_PROTECT_LAYER;

%First ebeam: ((NegLayer - SiLayer) OR HoleLayer)


%% GLOBAL PARAMETERS
% P for Parameters
P.write_field_x = 60;
P.write_field_y = 60;


P.bragg_point_res = 0.010; % use fine resolution for CAD only, not for plotting


P.trench = 3; % trench around normal waveguides
P.taper_length = 345; % for grating couplers
P.total_wg_length = 270 + 290; % length between end of GC tapers (default: 400)
P.scaling_pitch = 1.0;  %dev16 from v1 chip
P.tag_offset = [-10;-100];
P.wg_horz_spacing = 620; %P.write_field_x;  %space waveguides by the write field
P.fishbone_wg_width = 0;%0.372;
P.ring_gc_horz_spacing = 46;
P.wg_width = 0.476; % corrected, should come out around 450 nm.
P.N_points = 3;  %this is for the RingRes function. NEED TO REMOVE FROM FNC!
P.contact_wire_width = 1; % width of wire protruding from slot mode converter
P.contact_wire_length = 80;
P.GC_hole_corr = false;
P.si_protect_width = 1.5;% half the trench is optimal
P.scaling_pitch = 1.0;  % (-) Scaling of the pitch
P.scaling_holes = 1.0; % (-) Scaling of the grating holes on top of potential correction factor
%EOS6_elec_wire_length = 200;
%Bandage_width=40;
%Bandage_height=202.5;

N_fins = 12;
P.fin_width = [0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.110,0.110];
P.fin_gap_fishbone=0.00; %distance fin to waveguide
P.fin_gap_electrode=2.0; %distance fin to left electrode
P.fin_length=[500.0,500.0,500,500.0,500.0,500,500.0,500.0,500,500.0,500.0,500]; %device.cavity_length - space for tapers
P.taper_length = 9; %quarter period of cosine
P.taper_amplitude = 41.6; %needs to be adjusted if change fin_width or taper_length
P.fin_BS_gap = [0.1,0.1,0.1,0.1,0,0.1,0.1,0.1,0,0.1,0.1,0.1];
%P.BS_Lint=23.07*exp(4.3*P.fin_BS_gap);
P.BS_Lint=[35.5,40,45,49,0,52,20,63,0,71,80,90];
P.cav_tunnel=[1,3,5,7,8,2,4,6,9,7.5,4.5,8.5];
P.cav_tunnel2=[2.5,1.5,3.5,0.1,5,8,3,2,1.5,7,0.5,4];

P.ring_radius=25.0;
P.fin_ring_gap=[0.100,0.110,0.120,0.090,0.105,0.105,0.115,0.085,0.100,0.105,0.110,0.100]; %top
P.fin_ring_gap_bot=[0.100,0.110,0.120,0.090,0.115,0.105,0.115,0.085,0.095,0.105,0.110,0.100]; %bottom same as top except fotr devices without BS indexes (5 and 9)
P.ring_width=0.476;
P.ring_height=50;
P.wg_ring_distance=0.3;

P.number_electrodes=[400,40,400,400,400,400,400,400,400,400,400,400];
P.electrode_pattern_period=[1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0];
P.electrode_tot_length=P.number_electrodes.*P.electrode_pattern_period; 
P.electrode_pos_y = P.fin_length./2-P.electrode_tot_length/2;%-EOS6_elec_wire_length; %how much of the fin sticks out at the bottom of the electrodes
P.fin_ground_gap= [0.150,0.300,0.200,0.110,0.150,0.200,0.300,0.250,0.200,0.140,0.150,0.150];
P.electrode_pattern_fin_gap=[0.080,0.150,0.100,0.075,0.100,0.120,0.200,0.050,0.125,0.090,0.040,0.070];
P.electrode_pattern_width=0.350;
%% Initialize
master_group    =  Group([0,0],[1;0],[0;1],{});

% pattern must be away from [0,0] to avoid errors in AutoCAD, so offset
% origin
origin = [0;0];

curr_x = origin(1) + P.wg_horz_spacing/2;
curr_y = origin(2);

device_counter = 1;
fishbone_locations = {};
fishbone_total_length = []; 


%% make the fishbone array


fishbone_grp = Group([0,0],[1;0],[0;1],{});

% want to shift center wavelength from 1540 to 1578, so new period is 434
% nm


%DEVICE PARAMETERS
wg_width = P.fishbone_wg_width;
cavity_length = 450*ones(1, N_fins);
% reflector parameters
N_ramp_in = 10 * ones(1, N_fins);
N_mirror = ones(1,N_fins);
N_ramp_out = 10 * ones(1, N_fins);
ampl = 0.100 * ones(1, N_fins);
period = 0.436;

device = {};
curr_x = curr_x - P.wg_horz_spacing/2;
for ii = 1:N_fins

    curr_x = curr_x + P.wg_horz_spacing;
    if ii==2
       %curr_x=curr_x-P.fin_ground_gap-P.fishbone_wg_width/2+P.fin_gap_electrode(2)+P.fin_width+P.fin_gap_fishbone;
       curr_x=curr_x-P.fin_ground_gap(ii)-P.fishbone_wg_width/2;
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



    fishbone_locations{end+1} = [curr_x;curr_y]; % keep track of locations for aligning electrodes.
    fishbone_total_length(end+1) = device.cavity_length + 2*device.period*(device.N_ramp_in + device.N_mirror + device.N_ramp_out);

    device_counter = device_counter + 1;
    if ii==2
       %curr_x=curr_x+P.fin_ground_gap+P.fishbone_wg_width/2-P.fin_gap_electrode(2)-P.fin_width-P.fin_gap_fishbone;
       curr_x=curr_x+P.fin_ground_gap(ii)+P.fishbone_wg_width/2;
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
cavity_grp = Group([0,0],[1;0],[0;1],{});
label_grp=Group([0,0],[1;0],[0;1],{});
%fin_grp.layer=L.SI_LAYER;

for ii=1:N_fins
    %Fin:
    
    curr_x = fishbone_locations2(1,ii);
    curr_y = fishbone_locations2(2,ii);

    fin=Rect(curr_x-P.fin_width(ii),curr_y-P.fin_length(ii)/2,P.fin_width(ii),P.fin_length(ii),L.SI_LAYER); %position x, position y, width, height
    fin_grp.addelement(fin);
    fin_neg=Rect(curr_x-P.fin_width(ii)-P.trench,curr_y-P.fin_length(ii)/2-P.trench,P.fin_width(ii)+2*P.trench,P.fin_length(ii)+2*P.trench,L.NEG_LAYER);
    fin_grp.addelement(fin_neg);
    fin_si_protect=Rect(curr_x-P.fin_width(ii)-(P.trench-P.si_protect_width),curr_y-P.fin_length(ii)/2-(P.trench-P.si_protect_width),P.fin_width(ii)+2*(P.trench-P.si_protect_width),P.fin_length(ii)+2*(P.trench-P.si_protect_width),L.HB_SILICON_PROTECT_LAYER);
    fin_grp.addelement(fin_si_protect);
    
    
    f2c.wg_width = P.fin_width(ii);
    %f2c.c1_radius=20.0;
    f2c.c2_radius=20.0;
    f2c.c3_radius=20.0;
    %f2c.c4_radius=20.0;
    %f2c.w1_length=100;
    f2c.w2_length=0.5;
    f2c.w3_length=20;
    f2c.w4_length=0.5;
    %f2c.w5_length=50;
    fin2cav=makeFin2FinCavity( curr_x,curr_y, f2c, P, L, B);
    fin2cav.mirror([curr_x;curr_y],[0:1]);
    fin2cav.translate([- f2c.w3_length - f2c.c2_radius - f2c.c3_radius + f2c.wg_width/2;curr_y + P.fin_length(ii)/2 + f2c.c2_radius + f2c.c3_radius + f2c.w2_length/2 + f2c.w4_length - f2c.wg_width]);
    fin_grp.addelement(fin2cav);
    
    fin2cav2=makeFin2FinCavity( curr_x,curr_y, f2c, P, L, B);
    fin2cav2.mirror([curr_x;curr_y],[0:1]);
    fin2cav2.translate([- f2c.w3_length - f2c.c2_radius - f2c.c3_radius + f2c.wg_width/2;curr_y + P.fin_length(ii)/2 + f2c.c2_radius + f2c.c3_radius + f2c.w2_length/2 + f2c.w4_length - f2c.wg_width]);
    %fin_grp.addelement(fin2cav);
    
    
    %Make cavity
    [sillyr,dellyr,grtlyr,hlelyr,crvlyr,sltlyr] = nsymair(P.cav_tunnel(ii));
    cav_grp = Group([0,0],[1;0],[0;1],{});
    cav_grp.addelement(sillyr);
    cav_grp.addelement(dellyr);
    cav_grp.addelement(grtlyr);
    cav_grp.addelement(hlelyr);
    cav_grp.addelement(crvlyr);
    cav_grp.addelement(sltlyr);
    cav_grp.rotate(90);
    cav_grp.translate([curr_x - 91.439;curr_y - 196.736 - 2.25  + P.fin_length(ii)/2 + f2c.c2_radius + f2c.c3_radius + f2c.w2_length/2 + f2c.w4_length - f2c.wg_width]);
    cavity_grp.addelement(cav_grp);
    
    
    [sillyr,dellyr,grtlyr,hlelyr,crvlyr,sltlyr] = nsymair(P.cav_tunnel2(ii));
    cav_grp2 = Group([0,0],[1;0],[0;1],{});
    cav_grp2.addelement(sillyr);
    cav_grp2.addelement(dellyr);
    cav_grp2.addelement(grtlyr);
    cav_grp2.addelement(hlelyr);
    cav_grp2.addelement(crvlyr);
    cav_grp2.addelement(sltlyr);
    cav_grp2.rotate(90);
    cav_grp2.translate([curr_x - 91.439;curr_y - 196.736 - 2.25  + P.fin_length(ii)/2 + f2c.c2_radius + f2c.c3_radius + f2c.w2_length/2 + f2c.w4_length - f2c.wg_width]);
    

    
    bot_cav_grp = Group([0,0],[1;0],[0;1],{});
    bot_cav_grp.addelement(cav_grp2);
    bot_cav_grp.addelement(fin2cav2);
    bot_cav_grp.rotate([curr_x - 91.439;curr_y - 196.736 - 2.25  + P.fin_length(ii)/2 + f2c.c2_radius + f2c.c3_radius + f2c.w2_length/2 + f2c.w4_length - f2c.wg_width],180);
    %bot_cav_grp.translate([curr_x;curr_y])
    bot_cav_grp.translate([P.fin_width(ii)/2 + 172 + 10.713;curr_y - P.fin_length(ii)/2 + 73 - 6.2 - 0.108]);

    cavity_grp.addelement(bot_cav_grp);
    
%     device.tag_text=num2str(ii);
%     if B.include_label
%       tag_pos = [curr_x; curr_y  - P.total_wg_length/2 - P.taper_length ] + P.tag_offset ;
%       label_grp.addelement(makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER));
%     end
  
    
    

%     %Make Directional coupler
%     if false && (ii==1 || ii==2 || ii==3 || ii==4 || ii==6 || ii==7 || ii==8 || ii==10 || ii==11 || ii==12)
% 
%         BS.wg_width = P.fin_width(ii);
%         BS.c1_radius=20.0;
%         BS.c2_radius=20.0;
%         BS.c3_radius=20.0;
%         BS.c4_radius=20.0;
%         BS.w1_length=100;
%         BS.w2_length=P.BS_Lint(ii);
%         BS.w3_length=10;
%         BS.w4_length=750;
%         BS.w5_length=50;
%         fin2=makeMechDirectionalCoupler( curr_x, curr_y, BS, P, L, B);
%         fin2.mirror([curr_x;curr_y],[0;1]);
%         fin2.translate([P.fin_BS_gap(ii) + P.fin_width(ii)/2;curr_y - P.fin_length(ii)/2 + P.electrode_pos_y(ii) - P.ring_height - P.BS_Lint(ii)/2]);
%         fin_grp.addelement(fin2);
%     end



end

master_group.addelement(cavity_grp);
master_group.addelement(fin_grp);




    


%% Make chip dicing marks


dice_mark_grp = Group([0,0],[1;0],[0;1],{});

diced_chip_width = 10000;
diced_chip_height = 5000;

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

rail_offset = 1200; 


P.rail_width = 700;
P.metal_wire_width = 10;
P.alignment_slop = 15;


pad_horz_offset = P.inner_maxX - 300 - centerX;


if B.split_rails
    % make rails
    rail_gap = 150;
    TR = [centerX - rail_gap/2; -rail_offset];
    BL = [centerX - pad_horz_offset ; -rail_offset - P.rail_width];%+ P.big_pad_width

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

    TR = [centerX + pad_horz_offset ; -rail_offset];%- P.big_pad_width
    BL = [centerX - pad_horz_offset ; -rail_offset - P.rail_width];%+ P.big_pad_width
    top_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    top_rail.mirror([centerX; centerY], [1;0])
    electrode_wires_grp.add( top_rail);
    
    pad_horz_offset = P.inner_maxX - 300 - 1400 - centerX;
    P.rail_width2=100;
    TR = [centerX + pad_horz_offset ; -rail_offset];%- P.big_pad_width
    BL = [centerX - pad_horz_offset ; -rail_offset - P.rail_width2];%+ P.big_pad_width

    bottom_rail = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    electrode_wires_grp.add( bottom_rail);
    
    TR = [centerX - pad_horz_offset ; -rail_offset];
    BL = [centerX - pad_horz_offset - 1400 ; -rail_offset - P.rail_width];
    left_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    electrode_wires_grp.add( left_pad);
    
    TR = [centerX + pad_horz_offset + 1400; -rail_offset];
    BL = [centerX + pad_horz_offset ; -rail_offset - P.rail_width];
    right_pad = CornerRect(BL,TR, L.MODULATION_ELEC_EOS3_LAYER);
    electrode_wires_grp.add( right_pad);
end

% make modulation electrodes.

for ii = 1:length(fishbone_locations)

    curr = fishbone_locations{ii};
    curr_x = curr(1);
    curr_y = curr(2);

    % % make EOS6 electrode part
    EOS6_elec_wire_width = 5;
    EOS6_elec_wire_length = 100;

    % define some points (to simplify making electrode structure) (making
    % right side of electrode)
    Ax1 = P.fishbone_wg_width/2 + P.fin_ground_gap(ii);

    Ay1 = P.electrode_tot_length(ii)/2;
    Ay1_arm = 200;
    Ax2 = Ax1 + EOS6_elec_wire_width;
    Ay2 = EOS6_elec_wire_width/2;
    Ax3 = EOS6_elec_wire_length - EOS6_elec_wire_width/2;
    Ax4 = EOS6_elec_wire_length + EOS6_elec_wire_width/2;

    points = {[Ax1;Ay1], [Ax2;Ay1], [Ax2;Ay2], [Ax3;Ay2], [Ax3;Ay1_arm], [Ax4;Ay1_arm], [Ax4;-Ay2], [Ax2;-Ay2], [Ax2;-Ay1], [Ax1;-Ay1]};

    TR_EOS6_electrode = Polygon(points, L.MODULATION_ELEC_EOS6_LAYER);
    TR_EOS6_electrode.translate([curr_x + P.fin_gap_electrode - P.fin_ground_gap(ii) ;curr_y-P.fin_length(ii)/2+P.electrode_tot_length(ii)/2+P.electrode_pos_y(ii)]);
    electrode_wires_grp.add(TR_EOS6_electrode);
    BR_EOS6_electrode = Polygon(points, L.MODULATION_ELEC_EOS6_LAYER);
    BR_EOS6_electrode.translate([curr_x-P.fin_ground_gap(ii)+P.fin_gap_electrode+P.fin_width(ii)+P.fin_gap_fishbone;curr_y+P.fin_length(ii)/2-P.electrode_tot_length(ii)/2-P.electrode_pos_y(ii)]);

    BR_EOS6_electrode.rotate([curr_x;curr_y], 180);
    electrode_wires_grp.add(BR_EOS6_electrode);

    %Patterning
    %left
    TC_pattern=[P.fin_gap_electrode - P.electrode_pattern_fin_gap(ii);P.electrode_pattern_width];
    BC_pattern=[0;0];
     
    for jj=1:P.number_electrodes(ii)
        finger=CornerRect(BC_pattern,TC_pattern, L.MODULATION_ELEC_EOS6_LAYER);
        finger.translate([curr_x-P.fin_width(ii)-P.fin_gap_electrode;curr_y+(jj-1)*P.electrode_pattern_period(ii)-P.fin_length(ii)/2+P.electrode_pos_y(ii)]);
        electrode_wires_grp.add(finger);
    end
    %right
    TC_pattern=[P.fin_gap_electrode - P.fin_ground_gap(ii);P.electrode_pattern_width];
    BC_pattern=[0;0];
     
    for jj=1:P.number_electrodes(ii)
        finger=CornerRect(BC_pattern,TC_pattern, L.MODULATION_ELEC_EOS6_LAYER);
        finger.translate([curr_x + P.fin_ground_gap(ii);curr_y+(jj-1)*P.electrode_pattern_period(ii)-P.fin_length(ii)/2+P.electrode_pos_y(ii)]);
        electrode_wires_grp.add(finger);
    end

    % % make EOS3 electrode part
    EOS3_elec_wire_width = 20;

    TRT = [Ax3 + EOS3_elec_wire_width/2; rail_offset+P.fin_length(ii)/2-P.electrode_tot_length(ii)/2-P.electrode_pos_y(ii)];
    BLT = [Ax3 - EOS3_elec_wire_width/2; -Ay2];

%     TRB = [Ax3 + EOS3_elec_wire_width/2; rail_offset+P.fin_length/2-P.electrode_tot_length/2-P.electrode_pos_y];
%     BLB = [Ax3 - EOS3_elec_wire_width/2; -Ay2+P.fin_length/2];
    TRB = [Ax3 + EOS3_elec_wire_width/2; Ay2-P.fin_length(ii)/2+P.electrode_tot_length(ii)/2+P.electrode_pos_y(ii)];
    BLB = [Ax3 - EOS3_elec_wire_width/2; -rail_offset];
    TR_EOS3_electrode = CornerRect(BLT,TRT, L.MODULATION_ELEC_EOS3_LAYER);
    TR_EOS3_electrode.translate([curr_x;curr_y-P.fin_length(ii)/2+P.electrode_tot_length(ii)/2+P.electrode_pos_y(ii)]);
    electrode_wires_grp.add(TR_EOS3_electrode);
    BL_EOS3_electrode = CornerRect(BLB,TRB, L.MODULATION_ELEC_EOS3_LAYER);
    BL_EOS3_electrode.translate([curr_x+2*Ay2-2*EOS6_elec_wire_length;curr_y]);
    %BL_EOS3_electrode.rotate([curr_x;curr_y], 180);
    electrode_wires_grp.add(BL_EOS3_electrode);


end


master_group.addelement(electrode_wires_grp);



%% Add alignment marks

if B.include_align_marks

    align_mark_grp = Group([0,0],[1;0],[0;1],{});

    align_inner_width = 0.1;
    align_outer_width = 5;
    align_inner_length = 4;
    align_arm_length = 100; % = total length of one arm (half of full cross length).

    centerX = (P.inner_maxX + P.inner_minX) / 2;
    centerY = origin(2);

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

    my_tree = Group([0;0], [1;0], [0;1], process_tree(boxX,boxY,X,Y));
    my_tree.translate([P.inner_minX + 1300; P.minY + 500]);
    logo_grp.add(my_tree);

    logo_scl_factor = 0.86; % make logos the same size
    my_logo = Group([0;0], [1;0], [0;1], process_logo(boxX*logo_scl_factor,boxY*logo_scl_factor,X,Y));
    my_logo.translate([P.inner_minX + 1900; P.minY + 500]);
    logo_grp.add(my_logo);

    master_group.add(logo_grp);


    % make chip tag
     tag_pos = [P.inner_minX + 2500; P.minY + 500] ;
     logo_grp.add(makeTag(tag_pos, chip_tag, B, L.LOGO_LAYER, 'box_size', 10));


end

%% Make rectangle bounding box for metal protection

master_group.add( Rect([centerX;centerY], 10000, 5000, L.HB_METAL_GAP_LAYER) );



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
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; centerY], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length - 900; centerY], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length - 900; P.inner_minY + align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; P.inner_maxY - align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; P.inner_minY + align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length - 900; P.inner_maxY - align_arm_length], L.MARK_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

% add two layer align mark structures to determine alignment between
% EOS6 metal and EOS6 photonics write.
align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX - 2*align_arm_length; P.inner_maxY - align_arm_length], L.MARK_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );
align_mark_grp.addelement( TwoLayerNestedAlignmentMark([centerX + 2*align_arm_length; P.inner_minY + align_arm_length], L.MARK_LAYER, L.EOS6_METAL_LAYER, 0.5, 0.5, align_inner_length, 2*align_inner_length) );


% make a copy of all these alignment marks, but protect them from Si etch
align_inner_width = 20;
align_outer_width = 20;
align_inner_length = 4;
align_arm_length = 100;

align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % P
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length + 900; centerY], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % Q
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % R
align_mark_grp.addelement( NestedAlignmentMark([centerX; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % S
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length - 900; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BR
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_minX + align_arm_length + 900; P.inner_minY + align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % BL
align_mark_grp.addelement( NestedAlignmentMark([P.inner_maxX - align_arm_length - 900; P.inner_maxY - align_arm_length], L.HB_SILICON_PROTECT_LAYER, align_inner_width, align_outer_width, align_inner_length, align_arm_length) ); % TR

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
%%

function thru_wg2 = makeCurvedStripWaveguideDevice( curr_x, curr_y, GC_fn, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

%   if B.heidelberg_si_protect == true
%       width = 2*P.trench + device.wg_width - 2*P.si_protect_width;
%       elements{end+1} = Rect([0; 0], width, P.total_wg_length*0.7, L.HB_SILICON_PROTECT_LAYER);
%   end
  

  %g1
  elements{end+1} = GC_fn();
  elements{end}.rotate([curr_x; curr_y], 180);
  elements{end}.translate([curr_x-device.w1_length - device.c1_radius + device.wg_width/2; curr_y + (device.w2_length)/2 + device.c1_radius - device.wg_width/2]);
  
  %w1
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w1_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
  
  %c1
  elements{end+1}=makeArc(pi/2,device.c1_radius,device.wg_width,200,P.trench,L.SI_LAYER,L.NEG_LAYER);
  elements{end}.translate([curr_x - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2]);
  
  %w2
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w2_length, P.trench, L.SI_LAYER, L.NEG_LAYER);

  %c2
  elements{end+1}=makeArc(pi/2,device.c2_radius,device.wg_width,200,P.trench,L.SI_LAYER,L.NEG_LAYER);
  %elements{end+1}=Polygon(points,L.SI_LAYER);
  elements{end}.rotate([curr_x; curr_y], -90);
  elements{end}.translate([curr_x - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 + device.c2_radius - device.wg_width/2 - device.c3_radius + device.wg_width/2]);
  
  %w3
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w3_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2]);

  %c3
  elements{end+1}=makeArc(pi,device.c3_radius,device.wg_width,400,P.trench,L.SI_LAYER,L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length - device.c3_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c3_radius + device.wg_width/2]);
  
  %w4
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w4_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) + device.w4_length/2 + device.wg_width/2; curr_y - device.w2_length/2 - 2*device.c3_radius + device.wg_width - device.c3_radius + device.wg_width/2]);
  
  %g2
  elements{end+1} = GC_fn();
  elements{end}.translate([curr_x - (device.w3_length - device.w4_length) - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - 2*device.c3_radius + device.wg_width - device.c3_radius + device.wg_width/2]);
  
  
  if B.include_label
      tag_pos = [curr_x-device.w1_length - device.c1_radius + device.wg_width/2; curr_y + (device.w2_length)/2 + device.c1_radius - device.wg_width/2] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, device.tag_text, B, L.NEG_LAYER);
      elements{end+1} = Rect(curr_x-device.w1_length - device.c1_radius + device.wg_width/2 + P.tag_offset(1),curr_y + (device.w2_length)/2 + device.c1_radius - device.wg_width/2 + P.tag_offset(2),40,15,L.HB);
  end
  % make a Group out of the elements
  thru_wg2 = Group([x0,y0],[1;0],[0;1],elements);
end

function thru_wg2 = makeCurvedStripWaveguideDevice2( curr_x, curr_y, GC_fn, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;

%   if B.heidelberg_si_protect == true
%       width = 2*P.trench + device.wg_width - 2*P.si_protect_width;
%       elements{end+1} = Rect([0; 0], width, P.total_wg_length*0.7, L.HB_SILICON_PROTECT_LAYER);
%   end
  

  %g1
  elements{end+1} = GC_fn();
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - (device.w1_length + device.c1_radius) - device.c4_radius + device.wg_width; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length]);
  
  %w5
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w5_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w5_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
  
  %c4
  elements{end+1}=makeArc(pi/2,device.c4_radius,device.wg_width,1000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 180);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius]);
  
  %w1
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w1_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w1_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
  %c1
  elements{end+1}=makeArc(pi/2,device.c1_radius,device.wg_width,1000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2]);
  
  %w2
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w2_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w2_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x; curr_y]);
  
  %c2
  elements{end+1}=makeArc(pi/2,device.c2_radius,device.wg_width,1000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], -90);
  elements{end}.translate([curr_x - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 ]);
  

  %w3
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w3_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w3_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);

  %c3
  elements{end+1}=makeArc(pi/2,device.c3_radius,device.wg_width,1000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2 - device.c3_radius + device.wg_width/2]);
  
  %w4
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w4_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w4_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);

  %g2
  elements{end+1} = GC_fn();
  elements{end}.rotate([curr_x; curr_y], -90);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length]);
  

  if B.include_label
     if device.tag_down 
      tag_pos = [curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, [device.tag_text 'B'], B, L.NEG_LAYER);
      elements{end}.mirror([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length]+P.tag_offset - [P.tag_offset(1);0],[1;0])
      elements{end}.mirror([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length]+P.tag_offset - [P.tag_offset(1);0],[0;1])
      elements{end+1} = Rect(tag_pos(1) - 18,tag_pos(2),40,15,L.HB_SILICON_PROTECT_LAYER);
     elseif ~device.tag_down & ~device.tag_upC
      tag_pos = [curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, [device.tag_text 'A'], B, L.NEG_LAYER);
      elements{end+1} = Rect(tag_pos(1) - 2 ,tag_pos(2) - 15,40,15,L.HB_SILICON_PROTECT_LAYER);
     else
      tag_pos = [curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length] + P.tag_offset ;
      elements{end+1} = makeTag(tag_pos, [device.tag_text 'C'], B, L.NEG_LAYER);
      elements{end+1} = Rect(tag_pos(1) -2,tag_pos(2) -15,40,15,L.HB_SILICON_PROTECT_LAYER);
    end
  end
  % make a Group out of the elements
  thru_wg2 = Group([x0,y0],[1;0],[0;1],elements);
end

function Arc_polygon = makeArc(angle,radius,width,N,trench,si_protect,layer_Si,layer_neg,layer_si_protect)
    elements={};
    s_i=(0:angle/N:angle);
    x_i=radius*cos(s_i);
    y_i=radius*sin(s_i);
    points=[x_i(1);y_i(1)];
    for ii=2:length(s_i)
        points=horzcat(points,[x_i(ii);y_i(ii)]);
    end

    x_i2=(radius-width)*cos(s_i);
    y_i2=(radius-width)*sin(s_i);
    points2=[x_i2(1);y_i2(1)];
    for ii=2:length(s_i)
        points2=horzcat(points2,[x_i2(ii);y_i2(ii)]);
    end
    points2=horzcat(points2,[x_i(end);y_i(end)]);
    points=horzcat(points,fliplr(points2));
    points=horzcat(points,[x_i(1);y_i(1)]);
    elements{end+1} = Polygon(points,layer_Si);
    %
    s_i=(0:angle/N:angle);
    x_i=(radius+trench)*cos(s_i);
    y_i=(radius+trench)*sin(s_i);
    points=[x_i(1);y_i(1)];
    for ii=2:length(s_i)
        points=horzcat(points,[x_i(ii);y_i(ii)]);
    end

    x_i2=(radius-trench-width)*cos(s_i);
    y_i2=(radius-trench-width)*sin(s_i);
    points2=[x_i2(1);y_i2(1)];
    for ii=2:length(s_i)
        points2=horzcat(points2,[x_i2(ii);y_i2(ii)]);
    end
    points2=horzcat(points2,[x_i(end);y_i(end)]);
    points=horzcat(points,fliplr(points2));
    points=horzcat(points,[x_i(1);y_i(1)]);
    elements{end+1} = Polygon(points,layer_neg);
    %
    s_i=(0:angle/N:angle);
    x_i=(radius+(trench-si_protect))*cos(s_i);
    y_i=(radius+(trench-si_protect))*sin(s_i);
    points=[x_i(1);y_i(1)];
    for ii=2:length(s_i)
        points=horzcat(points,[x_i(ii);y_i(ii)]);
    end

    x_i3=(radius-(trench-si_protect)-width)*cos(s_i);
    y_i3=(radius-(trench-si_protect)-width)*sin(s_i);
    points2=[x_i3(1);y_i3(1)];
    for ii=2:length(s_i)
        points2=horzcat(points2,[x_i3(ii);y_i3(ii)]);
    end
    points2=horzcat(points2,[x_i(end);y_i(end)]);
    points=horzcat(points,fliplr(points2));
    points=horzcat(points,[x_i(1);y_i(1)]);
    elements{end+1} = Polygon(points,layer_si_protect);
    
    Arc_polygon = Group([0,0],[1;0],[0;1],elements);
end


function thru_wg2 = makeMechDirectionalCoupler( curr_x, curr_y, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;
  

  %w5
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w5_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w5_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
  
  %c4
  elements{end+1}=makeArc(pi/2,device.c4_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 180);
  elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius]);
  
  %w1
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w1_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w1_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
  %c1
  elements{end+1}=makeArc(pi/2,device.c1_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2]);
  
  %w2
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w2_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w2_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x; curr_y]);
  
  %c2
  elements{end+1}=makeArc(pi/2,device.c2_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], -90);
  elements{end}.translate([curr_x - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 ]);
  

  %w3
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w3_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w3_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);

  %c3
  elements{end+1}=makeArc(pi/2,device.c3_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2 - device.c3_radius + device.wg_width/2]);
  
  %w4
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w4_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w4_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);

  % make a Group out of the elements
  thru_wg2 = Group([x0,y0],[1;0],[0;1],elements);
end


function thru_wg2 = makeFin2FinCavity( curr_x, curr_y, device, P, L, B)
  elements = {};
  x0 = curr_x;
  y0 = curr_y; % save coordinates for shifting later

  curr_x = 0;
  curr_y = 0;
  

%   %w5
%   elements{end+1} = makeStripWaveguide(device.wg_width, device.w5_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
%   elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
%   elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w5_length, L.HB_SILICON_PROTECT_LAYER);
%   elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width - device.c4_radius; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius + device.w5_length/2]);
%   
%   %c4
%   elements{end+1}=makeArc(pi/2,device.c4_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
%   elements{end}.rotate([curr_x; curr_y], 180);
%   elements{end}.translate([curr_x - device.w1_length - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width + device.c4_radius]);
%   
%   %w1
%   elements{end+1} = makeStripWaveguide(device.wg_width, device.w1_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
%   elements{end}.rotate([curr_x; curr_y], 90);
%   elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
%   elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w1_length, L.HB_SILICON_PROTECT_LAYER);
%   elements{end}.rotate([curr_x; curr_y], 90);
%   elements{end}.translate([curr_x - device.w1_length/2  - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2 + device.c1_radius - device.wg_width/2]);
%   %c1
%   elements{end+1}=makeArc(pi/2,device.c1_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
%   elements{end}.translate([curr_x - device.c1_radius + device.wg_width/2; curr_y + device.w2_length/2]);
  
  %w2
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w2_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w2_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x; curr_y]);
  
  %c2
  elements{end+1}=makeArc(pi/2,device.c2_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], -90);
  elements{end}.translate([curr_x - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 ]);
  

  %w3
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w3_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w3_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length/2 - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2]);

  %c3
  elements{end+1}=makeArc(pi/2,device.c3_radius,device.wg_width,2000,P.trench,P.si_protect_width,L.SI_LAYER,L.NEG_LAYER,L.HB_SILICON_PROTECT_LAYER);
  elements{end}.rotate([curr_x; curr_y], 90);
  elements{end}.translate([curr_x - device.w3_length - device.c2_radius + device.wg_width/2; curr_y - device.w2_length/2 - device.c2_radius + device.wg_width/2 - device.c3_radius + device.wg_width/2]);
  
  %w4
  elements{end+1} = makeStripWaveguide(device.wg_width, device.w4_length, P.trench, L.SI_LAYER, L.NEG_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);
  elements{end+1}=Rect([0; 0], device.wg_width + 2*P.trench - 2*P.si_protect_width,device.w4_length, L.HB_SILICON_PROTECT_LAYER);
  elements{end}.translate([curr_x - (device.w3_length + device.c2_radius) - device.c3_radius + device.wg_width; curr_y - device.w2_length/2 - device.c3_radius + device.wg_width/2 - device.c2_radius + device.wg_width/2 - device.w4_length/2]);

  % make a Group out of the elements
  thru_wg2 = Group([x0,y0],[1;0],[0;1],elements);
end




function fin_ring_res = make_fin_ring_res(L_int_l,L_int_r,radius,width,N,trench,layer_Si,layer_neg,layer_hole,layer_inner_pos)
    elements={};
    angle1=tan(0.5*L_int_r/radius);
    angle2=pi-tan(0.5*L_int_l/radius);
    s_i=(angle1:(angle2-angle1)/N:angle2);
    x_i=radius*cos(s_i);
    y_i=radius*sin(s_i);
    points=[x_i(1);-y_i(1)];
    for ii=1:length(s_i)
        points=horzcat(points,[x_i(ii);y_i(ii)]);
    end
    for ii=1:length(s_i)
        points=horzcat(points,[x_i(end-ii+1);-y_i(end-ii+1)]);
    end
    elements{end+1} = Polygon(points,layer_Si);
    
    angle1=tan(0.5*L_int_r/(radius-width));
    angle2=pi-tan(0.5*L_int_l/(radius-width));
    s_i=(angle1:(angle2-angle1)/N:angle2);
    x_i=(radius-width)*cos(s_i);
    y_i=(radius-width)*sin(s_i);
    points2=[x_i(1);-y_i(1)];
    for ii=1:length(s_i)
        points2=horzcat(points2,[x_i(ii);y_i(ii)]);
    end
    for ii=1:length(s_i)
        points2=horzcat(points2,[x_i(end-ii+1);-y_i(end-ii+1)]);
    end
    elements{end+1} = Polygon(points2,layer_hole);
    
    elements{end+1}=Circ([0;0], radius - width/2 - trench, layer_inner_pos);
    elements{end+1}=Circ([0;0], radius + width/2 + trench, layer_inner_pos);


%     x_i2=(radius-width)*cos(s_i);
%     y_i2=(radius-width)*sin(s_i);
%     points2=[x_i2(1);y_i2(1)];
%     for ii=2:length(s_i)
%         points2=horzcat(points2,[x_i2(ii);y_i2(ii)]);
%     end
%     points2=horzcat(points2,[x_i(end);y_i(end)]);
%     points=horzcat(points,fliplr(points2));
%     points=horzcat(points,[x_i(1);y_i(1)]);
%     elements{end+1} = Polygon(points,layer_Si);
    
%     s_i=(0:angle/N:angle);
%     x_i=(radius+trench)*cos(s_i);
%     y_i=(radius+trench)*sin(s_i);
%     points=[x_i(1);y_i(1)];
%     for ii=2:length(s_i)
%         points=horzcat(points,[x_i(ii);y_i(ii)]);
%     end
% 
%     x_i2=(radius-trench-width)*cos(s_i);
%     y_i2=(radius-trench-width)*sin(s_i);
%     points2=[x_i2(1);y_i2(1)];
%     for ii=2:length(s_i)
%         points2=horzcat(points2,[x_i2(ii);y_i2(ii)]);
%     end
%     points2=horzcat(points2,[x_i(end);y_i(end)]);
%     points=horzcat(points,fliplr(points2));
%     points=horzcat(points,[x_i(1);y_i(1)]);
%     elements{end+1} = Polygon(points,layer_neg);
    fin_ring_res = Group([0,0],[1;0],[0;1],elements);
end



function [metagrating] = gratingcoupler_airSOI220_27deg_twostepmeta_focusing_fromRVL(silicon_layer,hole_layer,neg_layer,si_protect_layer,boolean_corrections,boolean_hide_holes,scaling_pitch,scaling_holes,wg_width,silicon_protect_width,write_field,varargin)
% Two-step grating coupler, angle 27 degrees, fully etched 220 SOI, not
% focusing (Benedikovic-style)
% FYI, scaling_pitch scales the longitudinal period and hole size, keeping duty cycle and transverse parameters fixed
% scaling_holes scales the hole sizes, keeping longitudinal and transverse
% pitch constant
% correction factors: keep periods the same (less than a percent deviation
% - in error bar), make holes smaller so they come out correctly
offset = [0 0];
metagrating = Group(offset(1),offset(2),{});
holes = Group(offset(1),offset(2),{},hole_layer);
%write_windows = Group(offset(1),offset(2),{},window_layer);
p = inputParser;
default_focusing_length = 37.8; % (um) Focusing length of the grating
default_defocus = 5.0; % (um) Testing the defocus
addParameter(p,'focusing_length',default_focusing_length,@isnumeric);
addParameter(p,'defocus',default_defocus,@isnumeric);
parse(p,varargin{:});
focusing_length = p.Results.focusing_length; % (um) Approximate desired focusing length, will be adjusted later
defocus = p.Results.defocus; % (um)

%% Hard-coded parameters
% check params: main period (810 OK), main duty cycle (0.55 OK), main SGW period (450 nm OK), main SWG gap (171 nm OK), match period (760 nm OK), match duty
% cycle (0.55 OK), match SWG period (450 nm OK), match SWG gap (124 nm OK), corrections for all (OK)
lambda0 = 1.55; % (um) Center wavelength
theta = 27; % (deg.) Grating angle in air at 1550 nm with perfect fab
NPz = 22; % (-) Number of periods along propagation
LambdaSWG = 0.45; % (um) Period transverse / SubWavelengthGrating SWG
DCz = 0.55; % (-) Duty cycle along propagation
W = 19; % (um) Width transverse
NPswg = round(W/LambdaSWG); % (-) Number of transverse holes
trench = 3; % (um)
Lambdaz = scaling_pitch*[0.76*ones(1,4),0.81*ones(1,18)]; % (um) Period along propagation

if boolean_corrections == 1
    % corrections are for dose 320, EOS6_1000pA April 2017
    % corrections from google doc "RVL - 201704 - Overview of chips"
    corr_holes_longi = 0.945; % OK
    corr_holes_longi_match = 0.956; % OK
    corr_holes_trans = 0.914; % OK
    corr_holes_trans_match = 0.869; % OK
    %corr_wire_width = 1.04; % no longer needed, is now an input argument
    %of this function
else
    corr_holes_longi = 1;
    corr_holes_longi_match = 1;
    corr_holes_trans = 1;
    corr_holes_trans_match = 1;
end

%% Deduced parameters
neff = lambda0/Lambdaz(5) + sind(theta); % (-) Effective refractive index in main section
neff_match = lambda0/Lambdaz(1) + sind(theta); % (-) Effective refractive index in matching section
focusing_length_match = focusing_length - sum(Lambdaz(1:4));
q_min = round(focusing_length*(neff - sind(theta))/lambda0); % (-) Minimal grating line number deduced from focusing length
focusing_length = q_min*lambda0/(neff - sind(theta)); % (um) Exact focusing length
q_min_match = round(focusing_length_match*(neff_match - sind(theta))/lambda0); % (-) Minimal grating line number deduced from focusing length
focusing_length_match = q_min_match*lambda0/(neff_match - sind(theta)); % (um) Exact focusing length
grating_length = sum(Lambdaz); % (um) Length of the metasection
L = focusing_length_match + defocus + grating_length; % (um) Length of grating plus focusing section
delta = atand(W/(2*L)); % (deg.) Divergence angle from focus point, keep around 14 degrees
LY = (L + trench)*tand(delta); % (um) Extra y-distance needed for trench behind grating


% First the grid of holes
G = [scaling_holes*corr_holes_trans_match*0.124*ones(1,4),scaling_holes*corr_holes_trans*0.171*ones(1,NPswg)]; % (um) Hole size transverse, corrected
Gz = Lambdaz*(1-DCz)*scaling_holes; % (um) Hole size longitudinal
Gz(1:4) = Gz(1:4)*corr_holes_longi_match; % Correct longitudinal hole size in matching section
Gz(5:end) = Gz(5:end)*corr_holes_longi; % Correct longitudinal hole size in main section

if boolean_hide_holes == 0
    for kk = 1:NPz
        for ll = 1:NPswg
            if kk <= 4 % if kk = 1 you're at the first grating line
                q = (kk-1) + q_min_match; % (-) Grating line number
                xshift = sum(Lambdaz(1:kk)) - Lambdaz(1);
                yshift = ((ll-1) - ((NPswg - 1)/2))*LambdaSWG - G(kk)/2;
                phi = atand(yshift/(xshift + focusing_length_match)); % (deg.) Polar angle of the hole, origin in focus
                %if abs(phi) <= delta + 1 % in case you want to get rid of
                %additional holes on top and bottom
                r = q*lambda0/(neff_match - cosd(phi)*sind(theta)) + 0.065 + defocus; % (um) Distance the hole should have to focus point. *** the 65 nm averages the period period in the matching and main section
                x = r*cosd(phi);
                y = r*sind(phi);
                hole = Rect(x,y,Gz(kk),G(kk),'base','corner');
                hole.rotate([x;y],phi);
                holes.add(hole);
            elseif kk >= 5 % if kk = 5 you're at the first main grating line
                q = (kk-5) + q_min; % (-) Grating line number
                xshift = sum(Lambdaz(1:kk)) - Lambdaz(1);
                yshift = ((ll-1) - ((NPswg - 1)/2))*LambdaSWG - G(kk)/2;
                phi = atand(yshift/(xshift - sum(Lambdaz(1:4)) + focusing_length)); % (deg.) Polar angle of the hole, origin in focus
                %if abs(phi) <= delta + 1 % in case you want to get rid of
                %additional holes on top and bottom
                r = q*lambda0/(neff - cosd(phi)*sind(theta)) + defocus; % (um) Distance the hole should have to focus point
                x = r*cosd(phi);
                y = r*sind(phi);
                hole = Rect(x,y,Gz(kk),G(kk),'base','corner');
                hole.rotate([x;y],phi);
                holes.add(hole);
            end
            %end
        end
    end
end

% Second the frame - origin is in center of first column of holes
etched_around_pts{1} = [0;wg_width/2 + trench];
etched_around_pts{2} = [L + trench;wg_width/2 + trench + LY];
etched_around_pts{3} = [L + trench;-wg_width/2 - trench - LY];
etched_around_pts{4} = [0;-wg_width/2 - trench];
hole_around = Polygon(etched_around_pts,neg_layer);

% Third add silicon layer for inner frame
silicon_inner_pts{1} = [0;-wg_width/2];
silicon_inner_pts{2} = [L;-W/2];
silicon_inner_pts{3} = [L;W/2];
silicon_inner_pts{4} = [0;wg_width/2];
silicon_inner = Polygon(silicon_inner_pts,silicon_layer);

% Fourth add silicon protect layer
silicon_protect_pts{1} = [0;-wg_width/2 - silicon_protect_width];
silicon_protect_pts{2} = [L + silicon_protect_width;-W/2 - silicon_protect_width];
silicon_protect_pts{3} = [L + silicon_protect_width;W/2 + silicon_protect_width];
silicon_protect_pts{4} = [0;wg_width/2 + silicon_protect_width];
silicon_protect = Polygon(silicon_protect_pts,si_protect_layer);

% Fifth add write window
%write_window = Rect(0,-write_field/2,write_field,write_field,'base','corner');
%write_windows.add(write_window);


%% Gather components and make the metagrating
metagrating.addelement(holes);
metagrating.addelement(hole_around);
metagrating.addelement(silicon_inner);
metagrating.addelement(silicon_protect);
%metagrating.addelement(write_windows);

end