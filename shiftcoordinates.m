% function xc = shiftcoordinates(dim,n,ncoord,x)
%
% Purpose
% =======
% Shift coordinates: xc = xi - x;
%
% Input variables
% ===============
% dim    = dimension of the domain of discretization
% n      = total number of node sin the domain
% ncoord = vector of nodal coordinates 
%          1d: ncoord = [x1;x2;...;xn]
%          2d: ncoord = [x1 y1;x2 y2;...;xn yn]
%          3d: ncoord = [x1 y1 z1;x2 y2 z2;...;x3 y3 z3]                                 
% x      = sampling (evaluation) point
%          1d: x = x or [x]
%          2d: x = [x y]
%          3d: x = [x y z]
%
% Return
% ======
% 1d: xc = [x1-x; x2-x;...; xn-x]
% 2d: xc = [x1-x y1-y; x2-x y2-y;...; xn-x yn-y]
% 3d: xc = [x1-x y1-y z1-z; x2-x y2-y z2-z;...; xn-x yn-y zn-z]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xc=shiftcoordinates(dim,n,ncoord,x)

if (dim == 1)
  xc=ncoord - [x(1)*ones(n,1)];
elseif (dim == 2)
  xc=ncoord - [x(1)*ones(n,1) x(2)*ones(n,1)];
elseif (dim == 3)
  xc=ncoord - [x(1)*ones(n,1) x(2)*ones(n,1) x(3)*ones(n,1)];
else
  error('Fatal error! Dimension not yet coded.')
end
