clear

addpath(genpath('/home/tim/github_local/CAD'))
addpath(genpath('/media/user_data/jwitmer/CAD_working_folder/2017_09_26_QuarterWaveQEOM_v1/'))

start_time = cputime;

master_group    =  Group([0,0],[1;0],[0;1],{});


%% BOOLEANS
B.generate_CAD = false;
B.use_high_fishbone_resolution = true;

B.include_label = true;
B.include_align_marks = true;
B.invert_tags = false;
B.hide_GC_holes = false;
B.include_logos = true;
B.include_dicing_marks = false;
B.split_rails = false;

%% make subchip A

% short couplers, with poly

L_coupler = 88;
chip_tag = 'A';
sub_group = QuarterWaveQEOM_v1_subchip(B, L_coupler, chip_tag);
ext = sub_group.getExtent;
BL = ext{1};
TR = ext{3};
sub_group.translate( -BL + [-3500; 3050]);
master_group.add( sub_group );


%% make subchip B

% long couplers, with poly

L_coupler = 214;
chip_tag = 'B';
sub_group = QuarterWaveQEOM_v1_subchip(B, L_coupler, chip_tag);
sub_group.translate( -BL + [-3500; 50]);
master_group.add( sub_group );


%% make subchip C

% short couplers, no silicon

L_coupler = 88;
chip_tag = 'C';
sub_group = QuarterWaveQEOM_v1_subchip(B, L_coupler, chip_tag);
sub_group.translate( -TR + [3500; -50]);
master_group.add( sub_group );

%% make subchip D

% air chip for silicon-only tests

electrode_gap = 1.2;
chip_tag = 'D';
sub_group = EOFishbone_v3_subchip(B, electrode_gap, chip_tag);
sub_group.translate([0;-4500] );
master_group.add( sub_group );


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
filename = 'QuarterWaveQEOM_v1_no_DieLayer';
display('Writing layers')
master_group.todxf(filename);
fclose('all')

end_time = cputime;

display(['Total time is: '  num2str(round(end_time - start_time)) ' seconds.'])
