addpath(genpath('/Users/Felix/Documents/GitHub/CAD/'));

a=0.080;
L=100;
etch=2;
%gap=0.100;
%R_oxide=3300;


silicon = gpack.Group(0,0,{},'silicon');
holes = gpack.Group(0,0,{},'holes');
%device = gpack.Group(0,0,{silicon,glass},'device');

% pointx=R_oxide*cos(pi*(0:0.001:1));
% pointy=-R_oxide*sin(pi*(0:0.001:1));
% point=vertcat(pointx,pointy);
% 
% rod3=Polygon(point);
% 
% glass.addelement(rod3);

%wg1=Rect(-gap/2-a,0,a,L);
%wg2=Rect(gap/2,0,a,L);

wg1=Rect(-etch-a/2,0,-a/2,L);
wg2=Rect(a/2,0,etch+a/2,L);

silicon.addelement(wg1);
silicon.addelement(wg2);

silicon.todxf('/Users/Felix/Documents/Master Thesis/Fab/Electromechanical_fin/CAD/generation1/phononWG');