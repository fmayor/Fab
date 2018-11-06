function [sillyr,dellyr,grtlyr,hlelyr,crvlyr,sltlyr] = nsymair(rs)
fname = 'nsymair';

% This pattern is based on FIN08a:fgcsym.
% It's documented in the parent fab dir's LOG.
%
% CJS
% July 29, 2016

% Coupling waveguide parameters 
%
wgW = 0.35;
wgL = 40;

% Buffers -- etched region around each device
%
buffer = 2;
centralBuffer = 4; % width of etch region near device (4 um anticipating functionalization of the fins)

% Grating coupler parameters
%
%
% Vertical grating couplers
%
% vertical grating coupler design by Yanni in his first GA run
% written on FIN09b
% tested by Joey
% Centers were at 1570 nm after gluing.  I'll scale the grating by 0.987 to
% center it at 1550 nm.
gc_h20_scale = -0.013;

fPath = '/Volumes/user_data/ydahmani/20160507_SiLN_GC/multi_variable_ga/';
fname1 = [fPath 'ga1Results_TE/vec/vec_12652.mat'];
gcw = 8;
tpw = 6;
gcL = 15;

load(fname1);
% Generate x_starts and x_ends and write into tmp_grating.mat
point_list = generate_point_list_from_dc(vec);
generate_grating_for_lumerical(point_list);
load('tmp_grating.mat');
% change units from m to um
x_starts = x_starts/1e-6;
x_ends = x_ends/1e-6;
% convert x_starts to widths of teeth and gratings
widths = x_ends(2:21) - x_starts(2:21);
fShift = x_starts(2) - x_ends(1);
x_pos = x_starts+(-x_starts(2)+fShift)*(1+gc_h20_scale);
widths = widths*(1+gc_h20_scale);

%
% Air grating couplers
%
% measured FIN09a device had a center at 1535
% if unit cell simulations are correct cavity should be at 1433 in air
gc_air_scale = 1.01;
gc_air_L = 20;
pitch_air = .628 * gc_air_scale;

sillyr = Group([0,0],[1;0],[0;1],{},'SiLayer');
%sillyr = Layer({},'name','SiLayer');
dellyr = Group([0,0],[1;0],[0;1],{},'NegLayer');
%dellyr = Layer({},'name','sublayer');
% gratings
grtlyr = Group([0,0],[1;0],[0;1],{},'GrtLayer');
%grtlyr = Layer({},'name','gratlayer');
% hole layer
hlelyr = Group([0,0],[1;0],[0;1],{},'HoleLayer');
%hlelyr = Layer({},'name','holes');
% curve lyr
crvlyr = Group([0,0],[1;0],[0;1],{},'CurveLayer');
%crvlyr = Layer({},'name','curves');
% slot box
sltlyr = Group([0,0],[1;0],[0;1],{},'slotbox');
%sltlyr = Layer({},'name','slotbox');
% bounding boxes
bblyr = Layer({},'name','boundingBoxes');

% no hole variation
depth = 1;

% FIN08a:sym:dev39
% fundamental mode at 1525 nm for 5% scaling
% sweeping around fundamental at 1530 (5.3% scaling)
%
% Simulations of nearly symmetric cavities were ~0.7% lower frequency for a
% second fin with constant 110 nm width.
% I shifted scales down from the FIN09:sym pattern.
%scales = 1.041:0.005:1.056;
scales=1.0460;
Nmirs=-14;
curve_geometries=[80 30 110 30]*1e-3;


% Nmir for kappa_e
%
% Studies done for symmetric fin cavities in air for 80 + 30 fins.
%
% Nmir =-18 -> Q_e = omega/kappa_e = 10k
% +8 unit cells -> Q_e +10 dB
%
%Nmirs = [-22:4:-10];

% Curve geoms
% 3 fin combos
% curve_geometries = [80 30 80  30;
%                     80 30 110 30;
%                     60 30 110 30]*1e-3;

inx = 0;

%% Through waveguide with TE and TM couplers

wgWs = [0.24,0.350];
%%
if false
for ii = [1:length(wgWs)]  
for jj = 2:2  % gratings

  %%
  inx = inx + 1;
  
  % extent of grating coupler
  if jj ==1 
    gc_extent = x_pos(21)+widths(end);
  elseif jj==2
    gc_extent = gc_air_L*pitch_air;
  end

  %% Construct waveguide

    y = inx*-62.5/2;

    wgW = wgWs(ii);
    
    % middle wg
      % parameters of geometry
      lc = 45;
      csegls = [0,lc,lc,lc,lc,3,wgL,3,lc,lc,lc,lc];
      cls = cumsum(csegls);
      ws = 0.5*[gcw,6,3,1,wgW,wgW];
      % wg
        points = genpointswg(gc_extent+buffer,y,cls,[ws flip(ws)]);
        wg     = Polygon(points);
      % trench layer
        ws(end) = ws(end) + 2;
        points = genpointswg(gc_extent+buffer,y,cls,[ws flip(ws)] + buffer);
        wgdel    = Polygon(points);

  %% Gratings
        
    x_gcl = buffer + gc_extent;
    x_gcr = buffer+gc_extent+cls(end);
    
   
    if jj == 1
      gc_left = GratingCoupler_f('x',x_gcl, 'y', y+gcw - gcw/2,'width',gcw, 'xPos', x_pos(2:21), 'widths', widths);
      gc_right = GratingCoupler_f('x',x_gcr,'y', y-gcw/2,'width',gcw, 'xPos', x_pos(2:21), 'widths', widths);
    elseif jj == 2
            
      duty_cycle_air = .85;
      firstTooth_scale = 1.06;
      
      % grating positions
      xpos_air = pitch_air*(1-duty_cycle_air):pitch_air:pitch_air*(gc_air_L-.85);
      x_init = pitch_air*(1-firstTooth_scale*duty_cycle_air);
      xpos_air(1) = x_init;
      % grating widths
      widths_air = pitch_air*duty_cycle_air*ones(1,gc_air_L);
      widths_air(1) = pitch_air*duty_cycle_air*firstTooth_scale;

      gc_left  = GratingCoupler_f('x',x_gcl,'y', y + gcw/2,'width',gcw, 'xPos', xpos_air, 'widths', widths_air);
      gc_right = GratingCoupler_f('x',x_gcr,'y', y - gcw/2,'width',gcw, 'xPos', xpos_air, 'widths', widths_air);
    end
      
    gc_left.rotate(180);
    gc_right.translate([0;4]);
    gc_left.translate([30-0.628;-1590+0.25]);
    
    % Layer from which gc_left and gc_right are substracted
    grtL = Rect(x_gcl - gc_extent,y-gcw/2,gc_extent,gcw);
    grtR = Rect(x_gcr,y-gcw/2,gc_extent,gcw);
    %grtR.translate([0;4]);
    
    % dellyr for gratings
    gcldel = Rect(x_gcl-buffer-gc_extent,-buffer+y - gcw/2,gc_extent + buffer,...
      gcw+2*buffer);
    gcrdel = Rect(x_gcr,-buffer + y - gcw/2,gc_extent*(1+gc_h20_scale) + buffer,gcw+2*buffer);


    total_length = gc_extent*2 + cls(end);

    % needed to avoid namespace collision later on  -- :/
    yoffset = y;
    clear y
    
  %% Add elements to layer

  sillyr.addelement(gc_left);
  sillyr.addelement(wg);
  sillyr.addelement(gc_right);

  dellyr.addelement(gcldel);
  dellyr.addelement(wgdel);
  dellyr.addelement(gcrdel);

  grtlyr.addelement(grtL);
  grtlyr.addelement(grtR);

   %% Number

%   if mod(inx,5) == 0
%     number = photonics_number(num2str(inx),'center',0);
%     number.translate(total_length+5,yoffset);
%     dellyr.addelement(number)
%   end


end
end
end
%%
%% Cavities
% n=1;
% for kk = [1:length(scales)]
% for ii = [1:size(curve_geometries,1)]
% for jj = [1:length(Nmirs)]
% disp(n)
% n=n+1;
% disp(scales(kk))
% disp(curve_geometries(ii,:))
% disp(Nmirs(jj))
% end
% end
% end



for kk = [1:length(scales)]
for ii = [1:size(curve_geometries,1)]
for jj = [1:length(Nmirs)]  
for gc_inx = 2:2

  %%
  inx = inx + 1;

  % extent of grating coupler
  if gc_inx ==1 
    gc_extent = x_pos(21)+widths(end);
  elseif gc_inx==2
    gc_extent = gc_air_L*pitch_air;
  end

  %% Make cavities

  %% Magic factors for dose 260/320
  
  mf_rx = 1.0865;
  mf_ry = 0.8955;
  mf_w  = 1.0572;
  mf_g  = 0.8421;
    
  %% Curve geometry
  
  w01 = curve_geometries(ii,1)+curve_geometries(ii,2);
  curve_d = curve_geometries(ii,2);
  curve_w = 15; %length of curve
  
  w02 = curve_geometries(ii,3)+curve_geometries(ii,4);
  curve_d2 = curve_geometries(ii,4);
  curve_w2 = 15;
  
  %% Holes
  
  load /Volumes/user_data/sicamor/numerics/soeye/TEcavity/symmetrizedUC/singlemoded/ga3_optcell/optcellgeom.mat

  if Nmirs(jj) < 0
    Nt = 40 + Nmirs(jj);
    Nm = 0;
  else
    Nt = 40;
    Nm = Nmirs(jj);
  end
  Nramp = 10;

  % scalings from xpoint sweep of a in this directory
  s_def = 1.00;
  s_mir = depth;
  
  sii = scales(kk);

  % since rx = ry for this defect cell and mirror cell, we can just set
  % theta to 0.
  fin_w = max(w01,w02);
  theta  = 0;
  defect = struct('a',a*s_def*1e-3*sii,'w',w*1e-3*mf_w*sii,...
                  'e_a',max(rx,ry)*1e-3*mf_rx*sii,'e_b',min(rx,ry)*1e-3*mf_ry*sii,...
                  'theta',theta,'e_x',0,'e_y',0,'g',60e-3*mf_g,'fin_w',fin_w);

  theta  = 0;
  mirror = struct('a',a*s_mir*1e-3*sii,'w',w*1e-3*mf_w*sii,...
                  'e_a',max(ry*s_mir,rx*s_mir)*1e-3*mf_rx*sii,'e_b',min(ry*s_mir,rx*s_mir)*1e-3*mf_ry*sii,...
                  'theta',theta,'e_x',0,'e_y',0,'g',60e-3*mf_g,'fin_w',fin_w);
  
  params = generate_holeParams(defect,mirror,Nt,Nm);
  params = [flip(params) params(2:end)];
  % add adiabatic ramp up and down
  rampUp(Nramp) = mirror;
  thetas = linspace(0,1,Nramp+2);
  thetas = thetas(2:end-1);
  for mm = 1:Nramp
    rampUp(mm) = params(end);
    rampUp(mm).e_a = rampUp(mm).e_a*thetas(mm);
    rampUp(mm).e_b = rampUp(mm).e_b*thetas(mm);
  end
  params = [rampUp params flip(rampUp)];

  [slab,holes] = make_cavity(@(P)cad_symfin_uc(P,'theta',90),params,'rotbeam',1); 

   
  cav_exx = sum([params.a]);
    
  slotbox = Rect(0,-mirror.w/2-mirror.g-mirror.fin_w,cav_exx,mirror.w+2*mirror.g+2*mirror.fin_w);

  if cav_exx > wgL
    error('cavity is longer then center waveguide');
  end
  
  %% Make shapes to be deleted from the fins to curve them
    
  N_curve = 1500; % number of points defining the curved region
  
  % points for polygon that define curve of the first (top) fin
  x = [linspace(-curve_w/2,curve_w/2,N_curve)];
  y = [curve_d*(2*x/curve_w).^2, curve_d, curve_d + w01, curve_d + w01,curve_d];
  x = [x, cav_exx/2, cav_exx/2, -cav_exx/2, -cav_exx/2];
  % points for polygon that define curve of the second (bottom) fin
  x2 = [linspace(-curve_w2/2,curve_w2/2,N_curve)];
  y2 = [curve_d2*(2*x2/curve_w2).^2, curve_d2, curve_d2 + w02, curve_d2 + w02,curve_d2];
  x2 = [x2, cav_exx/2, cav_exx/2, - cav_exx/2, -cav_exx/2];
  %rs=5;
  s=cav_exx/2-rs;

  thick_diff=curve_geometries(4); %0.03
  x2 = [linspace(-curve_w2/2,curve_w2/2,N_curve)];
  y2 = [curve_d2*(2*x2/curve_w2).^2];
  y2 = [y2,curve_d2, curve_d2 + w02, curve_d2 + w02,curve_d2, y2(1)];
  x2 = [x2,cav_exx/2, cav_exx/2, - s, -s,x2(1)];
 
  
  y3 = [curve_d2, curve_d2 + w02, curve_d2 + w02,curve_d2,curve_d2];
  x3 = [- cav_exx/2 + rs, - cav_exx/2 + rs, - cav_exx/2, -cav_exx/2,- cav_exx/2 + rs];
  


  a4= (y2(end)-(y2(end-1)-thick_diff))/(x2(end)-x2(end-1))^2;
  b4=-2*a4*x2(end-1);
  c4=y2(end-1)-thick_diff + a4*x2(end-1)^2;
  
  xtmp = [linspace(x2(end-1),x2(end),N_curve)];
  x4 = [xtmp,xtmp(1),xtmp(1)];
  y4 = [a4*xtmp.^2+b4*xtmp+c4,y2(end-1),y2(end-1) - thick_diff];
  
  figure()
  plot(x2,y2,x3,y3-thick_diff,x4,y4)
  
  
  top_curve = Polygon(arrayfun(@(x,y)[x,y],x,y,'UniformOutput',false)); 
  top_curve.translate(cav_exx/2,mirror.w/2+mirror.g+w01 - curve_d);
  
  bot_curve = Polygon(arrayfun(@(x,y)[x,-y],x2,y2,'UniformOutput',false));
  bot_curve.translate(cav_exx/2,-(mirror.w/2+mirror.g+w02 - curve_d2));
  
  bot_curve2 = Polygon(arrayfun(@(x,y)[x,-y],x3,y3,'UniformOutput',false));
  bot_curve2.translate(cav_exx/2,thick_diff -(mirror.w/2+mirror.g+w02 - curve_d2));
  
  bot_curve3 = Polygon(arrayfun(@(x,y)[x,-y],x4,y4,'UniformOutput',false));
  bot_curve3.translate(cav_exx/2,-(mirror.w/2+mirror.g+w02 - curve_d2));
  
  curves = Group(0,0,{top_curve,bot_curve,bot_curve2,bot_curve3});

  %% Construct waveguide

    y = inx*-62.5/2;

    wgW = mirror.w;
    
    % middle wg
      % parameters of geometry
        lc = 45;
        csegls = [0,lc,lc,lc,lc,3,wgL,3,lc,lc,lc,lc];
        cls = cumsum(csegls);
        ws = 0.5*[gcw,6,3,1,0.35,wgW];
      % wg
        points = genpointswg(gc_extent+buffer,y,cls,[ws flip(ws)]);
        wg     = Polygon(points);
      % trench layer
        ws(end) = ws(end) + 2;
        points = genpointswg(gc_extent+buffer,y,cls,[ws flip(ws)] + buffer);
        wgdel    = Polygon(points);

  %% Gratings
        
    x_gcl = buffer + gc_extent;
    x_gcr = buffer+gc_extent+cls(end);
    
    if gc_inx == 1
      gc_left = GratingCoupler_f('x',x_gcl, 'y', y+gcw - gcw/2,'width',gcw, 'xPos', x_pos(2:21), 'widths', widths);
      gc_right = GratingCoupler_f('x',x_gcr,'y', y-gcw/2,'width',gcw, 'xPos', x_pos(2:21), 'widths', widths);
    elseif gc_inx == 2
            
      duty_cycle_air = .85;
      firstTooth_scale = 1.06;
      
      % grating positions
      xpos_air = pitch_air*(1-duty_cycle_air):pitch_air:pitch_air*(gc_air_L-.85);
      x_init = pitch_air*(1-firstTooth_scale*duty_cycle_air);
      xpos_air(1) = x_init;
      % grating widths
      widths_air = pitch_air*duty_cycle_air*ones(1,gc_air_L);
      widths_air(1) = pitch_air*duty_cycle_air*firstTooth_scale;

      gc_left  = GratingCoupler_f('x',x_gcl,'y', y + gcw/2,'width',gcw, 'xPos', xpos_air, 'widths', widths_air);
      gc_right = GratingCoupler_f('x',x_gcr,'y', y - gcw/2,'width',gcw, 'xPos', xpos_air, 'widths', widths_air);
    end
    gc_left.rotate(180);
    gc_right.translate([0;4]);
    gc_left.translate([30-0.628;-58-0.5]);


    % Layer from which gc_left and gc_right are substracted
    grtL = Rect(x_gcl - gc_extent,y-gcw/2,gc_extent,gcw);
    grtR = Rect(x_gcr,y-gcw/2,...
      gc_extent,gcw);
    
    % dellyr for gratings
    gcldel = Rect(x_gcl-buffer-gc_extent,-buffer+y - gcw/2,gc_extent + buffer,...
      gcw+2*buffer);
    gcrdel = Rect(x_gcr,-buffer + y - gcw/2,gc_extent*(1+gc_h20_scale) + buffer,gcw+2*buffer);


    total_length = gc_extent*2 + cls(end);

    % needed to avoid namespace collision later on  -- :/
    yoffset = y;
    clear y
    
    
  %% Add elements to layer for waveguide and gratings

  sillyr.addelement(gc_left);
  sillyr.addelement(wg);
  sillyr.addelement(gc_right);

  dellyr.addelement(gcldel);
  dellyr.addelement(wgdel);
  dellyr.addelement(gcrdel);

  grtlyr.addelement(grtL);
  grtlyr.addelement(grtR);

  %% Place the holes in the center of the waveguide.

  tranx = total_length/2 - cav_exx/2;
  trany = yoffset;

  holes.translate(tranx,trany);
  slab.translate(tranx,trany);
  curves.translate(tranx,trany);
  slotbox.translate(tranx,trany);
  hlelyr.addelement(holes);
  sillyr.addelement(slab);
  crvlyr.addelement(curves);
  sltlyr.addelement(slotbox);
  %% Number

%   if mod(inx,5) == 0
%     number = photonics_number(num2str(inx),'center',0);
%     number.translate(total_length+5,yoffset);
%     dellyr.addelement(number)
%   end


end
end
end
end

%% Collect Layers

% silicon layers
silyrs = {sillyr};
% trench layers
dellyrs = {dellyr,grtlyr,hlelyr,crvlyr,sltlyr};

lyrs = {silyrs{:}, dellyrs{:}};

%% Translate

%Translate everything upward and to the right
%foo_translate = @(y) cellfun(@(x) x.translate(0,15),y.elements);
%cellfun(foo_translate,lyrs);

%% Bounding Boxes

winl = 62.5;
brw = 9;
brh = ceil(inx/2)+2;
brx = -1;
bry = -(ceil(inx/2)+1);
boundingRectangle = Rect(brx*winl,bry*winl,brw*winl,brh*winl);  
bblyr.addelement(boundingRectangle.toAlignmentMarks());
lyrs{end+1} = bblyr;

%% Plot

% h_slotringres = figure;
% 
% clf
% drawnow
% 
% %%
% 
% set(0,'CurrentFigure',h_slotringres);
% for layerloop = silyrs
%   layerloop{1}.draw;
%   drawnow
% end
% 
% axis equal
% axis tight
% 
% drawwindowgrid();

%%


%write_file([fname ''],lyrs)

%save([fname '_data'],'boundingRectangle','fname');
end
