classdef GratingCoupler_f < Group
  %
  %
  % (x,y) position
  % length in number of unit cells
  % width in microns; transverse to propagation direction
  % pitch in microns; periodicity
  % dutyCycle : width of silicon / period
  %

  
  properties
    length
    pitchDutyArr
    widths
    xPos
    width
    pitch
    dutyCycle
    buffer
    direction
  end
  
  methods
    
    function [obj] = GratingCoupler_f(varargin)
      
      % parse arguments
      p = inputParser;
      padp = @(varargin) p.addParameter(varargin{:});
      padp('length',20);
      padp('width', 2);
      padp('dutyCycle',0.5);
      padp('pitchDutyArr',-1);
      padp('xPos',-1);
      padp('widths',-1);
      padp('pitch',0.5);
      padp('x',0);
      padp('y',0);
      padp('buffer',2);
      padp('direction','right');
      p.parse(varargin{:});
      P = p.Results;
      
      % instantiate group
      obj@Group(P.x,P.y,{});
      
      % initialize fields
      fnames = fieldnames(obj);
     % fieldnames(P)
      for i = 1:length(fnames)
        if ~strcmpi(fnames(i),'elements') && ~strcmpi(fnames(i),'x_axis') && ~strcmpi(fnames(i),'y_axis') && ~strcmpi(fnames(i),'origin') && ~strcmpi(fnames(i),'layer') && ~strcmpi(fnames(i),'is_hidden') && ~strcmpi(fnames(i),'name')
          obj.(fnames{i}) = P.(fnames{i});
        end
      end
      
      if obj.pitchDutyArr == -1 & obj.widths == -1
        for ii = 1:obj.length
            x = (ii-1)*obj.pitch;
            y = 0;
            w = obj.dutyCycle*obj.pitch;
            h = obj.width;
            obj.addelement(Rect(x,y,w,h, 'SiLayer', 0, 'base','corner'));
        end
      elseif obj.widths ~=-1
        for ii = 1:length(obj.widths)
            x = obj.xPos(ii);
            y = 0;
            w = obj.widths(ii);
            h = obj.width;
            obj.addelement(Rect(x,y,w,h,'SiLayer',0,'base','corner'));
        end
        obj.length = obj.xPos(end)+obj.widths(end);
      else
        pitchSum = [0 cumsum(obj.pitchDutyArr(1:2:end))];
        for ii = 1:2:length(obj.pitchDutyArr)
            jj = floor(ii/2)+1;
            x = (1-obj.pitchDutyArr(ii+1))*obj.pitchDutyArr(ii)+pitchSum(jj);
            y = 0;
            w = obj.pitchDutyArr(ii)*obj.pitchDutyArr(ii+1);
            h = obj.width;
            obj.addelement(Rect(x,y,w,h,  'SiLayer', 0, 'base', 'corner'));
        end
        obj.length = pitchSum(end);
      end
      if strcmpi(obj.direction,'right')
          obj.addelement(Rect(0,0,obj.length,obj.width,'GrtLayer',0,'base','corner'));
          obj.addelement(Rect(0, -obj.buffer,obj.length+obj.buffer,...
          obj.width+2*obj.buffer,'NegLayer',0,'base','corner'));
          obj.translate(0,-obj.width/2);
      else
          x = obj.x;
          y = obj.y;
          obj.translate(-x,-y);
          obj.rotate(180);
          obj.translate(x,y);
          %obj.translate(obj.length,obj.width);
          obj.addelement(Rect(-2*obj.length,-2*obj.width,obj.length,obj.width,'GrtLayer',0,'base','corner'));
      obj.addelement(Rect(-obj.buffer-2*obj.length,-2*obj.width -obj.buffer,obj.length+obj.buffer,...
          obj.width+2*obj.buffer,'NegLayer',0,'base','corner'));
      obj.translate(-obj.length,-obj.width/2);
          
      end
    end
    
  end
end