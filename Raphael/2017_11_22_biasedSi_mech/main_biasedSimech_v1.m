clear

%addpath(genpath('C:\Users\rvanlaer\Documents\GitHub\CAD\'))
%addpath(genpath('C:\Users\rvanlaer\Documents\GitHub\CAD\device\photonics\RVL_various_220SiPhot\'))
addpath(genpath('/Users/Felix/Documents/GitHub/CAD/'));
addpath(genpath('/Users/Felix/Documents/Master Thesis/Fab/Raphael/2017_11_22_biasedSi_mech/'));

start_time = cputime;

master_group    =  Group([0,0],[1;0],[0;1],{});


%% BOOLEANS
B.generate_CAD = true;
B.use_high_fishbone_resolution = false;

B.include_label = true;
B.include_align_marks = true;
B.invert_tags = false;
B.hide_GC_holes = false;
B.include_logos = true;
B.include_dicing_marks = false;
B.split_rails = false;

%% make subchip A

% short couplers, with poly

electrode_gap = 1.2;
chip_tag = 'A';
sub_group = EOFishbone_v3_subchip(B, electrode_gap, chip_tag);
% ext = sub_group.getExtent;
% BL = ext{1};
% TR = ext{3};
BL = 0;
TR = 0;
sub_group.translate( -BL + [-3500; 3050]);
master_group.add( sub_group );


%% make subchip B

% long couplers, with poly

% electrode_gap = 1.2;
% chip_tag = 'B';
% sub_group = EOFishbone_v3_subchip(B, electrode_gap, chip_tag);
% sub_group.translate( -BL + [-3500; 50]);
% master_group.add( sub_group );


%% make subchip C

% short couplers, no silicon

% electrode_gap = 1.2;
% chip_tag = 'C';
% sub_group = EOFishbone_v3_subchip(B, electrode_gap, chip_tag);
% sub_group.translate( -BL + [-3500; -2950]);
% master_group.add( sub_group );

%% make subchip D

% air chip for silicon-only tests

% electrode_gap = 1.2;
% chip_tag = 'D';
% sub_group = EOFishbone_v3_subchip(B, electrode_gap, chip_tag);
% sub_group.translate([-3500;-5950] );
% master_group.add( sub_group );


%% add rect for 10x15 mm die box

%master_group.add( Rect([0;0], 10000, 15000, 'DieLayer') );


%% Add extent marks: (every indepedent pattern should have same extents, makes aligned beamwrite much easier)

% NOTE: if a new layer is added in subchip file, must update this array to include all layers!
allLayers = {'SiLayer', 'NegLayer', 'HoleLayer', 'InnerPosLayer', ...
    'HBSiProtectLayer', 'HBMetalGapLayer', 'HBMetalPosLayer', 'EOS6MetalLayer', 'ExtentLayer'};

maxY = 6180 + 30; % want ring res location (y = -4500) to be at center of write field
minY = -6180 - 30;
minX = -3720 - 30;  % want x = 0 to be at center of write field
maxX = 3720 + 30;

extent_mark_size = 5;

for ii = 1:length(allLayers)
    %TL

    nextMark = CornerRect([minX;maxY - extent_mark_size], [minX + extent_mark_size; maxY], allLayers{ii});
    master_group.add(nextMark);

    %TR
    nextMark = CornerRect([minX;maxY - extent_mark_size], [minX + extent_mark_size; maxY], allLayers{ii});
    nextMark.mirror([0;0], [0;1]);
    master_group.add(nextMark);

    %BL
    nextMark = CornerRect([minX;maxY - extent_mark_size], [minX + extent_mark_size; maxY], allLayers{ii});
    nextMark.mirror([0;0], [1;0]);
    master_group.add(nextMark);

    %BR
    nextMark = CornerRect([minX;maxY - extent_mark_size], [minX + extent_mark_size; maxY], allLayers{ii});
    nextMark.mirror([0;0], [0;1]);
    nextMark.mirror([0;0], [1;0]);
    master_group.add(nextMark);


end



%% make DXF
% dxf file name, without .dxf
filename = '/Users/Felix/Documents/Master Thesis/Fab/Electromechanical_fin/CAD/generation1/Heidelwrite_RVL_goldtests_v_felix';
display('Writing layers')
master_group.todxf(filename);
fclose('all')

end_time = cputime;

display(['Total time is: '  num2str(round(end_time - start_time)) ' seconds.'])
