clear;

tic;

addpath(genpath('/home/tim/github_local/CAD'))
addpath(genpath('C:/Users/debussy/Documents/GitHub/CAD'))

%% BOOLEANS
B.include_label = true;
B.invert_tags = false;
B.hide_GC_holes = false;
B.heidelberg_si_protect = false;
B.generate_CAD = true;

if B.generate_CAD
    P.bragg_point_res = 0.010; % use fine resolution for CAD only, not for plotting
else
    P.bragg_point_res = 1.00;
end

%% PARAMETERS
P.write_field_x = 60;
P.write_field_y = 60;
P.trench = 3.0; % trench around normal waveguides
P.tag_offset = [-10;-45];
P.N_points = 3; %isnt used
P.font_size = 0.600;
P.taper_length = 250;
P.total_wg_length = 300; % length between end of GC tapers (default: 400)

P.wg_horz_spacing = P.write_field_x;  %space waveguides by the write field
P.wg_width = 0.300;

%% LAYERS
L.SI_LAYER = 'SiLayer';
L.NEG_LAYER = 'NegLayer';
L.HOLE_LAYER = 'HoleLayer'; % Si to be removed from positive regions.
L.INNER_POS_LAYER = 'InnerPosLayer'; % for ring inside
L.layer_names = {'SiLayer', 'NegLayer', 'HoleLayer', 'InnerPosLayer'};

first_group = Group([0,0],[1;0],[0;1],{});
second_group = Group([0,0],[1;0],[0;1],{});
origin = [0; 0];
curr_x = origin(1);
curr_y = origin(2);

dev_counter = 1;

%% do normal grating coupler sweep
[g, dev_counter] = normal_GC_sweep(P, L, B, dev_counter);
%first_group.add(g);

%% sweep fishbone cavities
[g, dev_counter] = fishbone_sweep(P, L, B, dev_counter);
%g.translate([1500; 0]);
%first_group.add(g);

%% sweep ring resonators
[g, dev_counter] = ring_sweep(P, L, B, dev_counter);
second_group.add(g);

%% add dot to set extent in the x-dimension
%dot = Rect(-30, 0, 1, 1,'base','corner','layer', L.NEG_LAYER);
%first_group.add(dot);

dot = Rect(-67, 390, 1, 1,'base','corner','layer', L.NEG_LAYER);
second_group.add(dot);

%first_group.todxf('PhotonRefresh_v1_1');
second_group.todxf('PhotonRefresh_v1_2');

toc
