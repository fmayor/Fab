function [slab, hole] = FBHZMaker_F(~, P)
% function [slab, hole] = FBHMaker(~, P)
% make fishbone-hole unitcell
% direction of periodicity is along x-axis
% WTJ, 20180712
% 
% Modified for FBH zipper nanobeam
% WTJ, FMM, 20181105
% 

N = 50;  % number of points for sin curve
w1 = P.w1+P.gap-0.150;
w2 = P.w2+P.gap-0.150;
w = (w1+w2)/2;
amp = P.amp;

if P.hx >= P.hy
   a = P.hx/2; 
   b = P.hy/2; 
   theta = 0;
else
    a = P.hy/2; 
    b = P.hx/2; 
    theta = 90; 
end
gap = P.gap;       % original gap, used for shifting the ellipse
gap_c = P.gap_c;   % corrected gap
%         slab = Rect(0,0,P.a,P.w,'base','center');
hole1 = Ellipse('rx',a,'ry',b,'theta',theta);
hole1.translate([0;gap/2]);
hole2 = Ellipse('rx',a,'ry',b,'theta',theta);
hole2.translate([0;-gap/2]);
w_slot = gap_c; h_slot = P.a;
rect_slot = Rect(0,0,h_slot,w_slot, 'base','center');
hole = gpack.Group(0,0,{hole1,hole2,rect_slot});

s = linspace(-pi, pi, N);
points = [s/2/pi * P.a; amp * cos(s) + w/2 + (w2-w1)/2 * s /2/pi];
points_flip = [flip(points(1,:)); flip(-points(2,:))];
points = [points points_flip];
slab = Polygon_F(points);

end