clear
addpath(genpath('/Users/Felix/Documents/localRepos/CAD'));
addpath(genpath('/Users/Felix/Documents/localRepos/Optomechanics'));

start_time = cputime;

master_group    =  Group([0,0],[1;0],[0;1],{});


%% BOOLEANS
B.generate_CAD = true;
B.include_label = true;
B.include_align_marks = true;
B.invert_tags = false;
B.hide_GC_holes = false;
B.include_logos = true;
B.include_dicing_marks = false;
B.split_rails = false;

%% make subchip A



chip_tag = 'A';
sub_group = Pilatus(B, chip_tag);
% ext = sub_group.getExtent;
% BL = ext{1};
% TR = ext{3};
BL = 0;
TR = 0;
%sub_group.translate( -BL + [-3500; 3050]);
master_group.add( sub_group );


%% add rect for 10x15 mm die box

%master_group.add( Rect([0;0], 10000, 5000, 'DieLayer') );


%% Add extent marks: (every indepedent pattern should have same extents, makes aligned beamwrite much easier)

% NOTE: if a new layer is added in subchip file, must update this array to include all layers!
allLayers = {'SiLayer', 'NegLayer', 'HoleLayer', 'InnerPosLayer', ...
    'HBSiProtectLayer', 'HBMetalGapLayer', 'HBMetalPosLayer', 'ExtentLayer','BandageLayer','SpiralLayer'};

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
filename = 'master';
display('Writing layers')
master_group.todxf(filename);
fclose('all')

end_time = cputime;

display(['Total time is: '  num2str(round(end_time - start_time)) ' seconds.'])
