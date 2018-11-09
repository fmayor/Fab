
function dev = genBondPad(x,y, w_bondPad)

if nargin < 3
    w_bondPad = 190;
end
    w_bondPadWindow = w_bondPad - 10;
    r_bondPad = Rect(0, 0, w_bondPad, w_bondPad, 'base','center');
    r_bondPad.layer = 'metal_BB';
    r_bondPadWin = Rect(0, 0, w_bondPadWindow, w_bondPadWindow, 'base','center');
    r_bondPadWin.layer = 'M2_undercutMsk';
    dev = gpack.Group(x,y,{r_bondPad,r_bondPadWin});
end
