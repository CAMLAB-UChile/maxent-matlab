% function h_node = nodespacing(dim,n,ncoord)
%
% Purpose
% =======
% Compute characteristic nodal spacing to use in basis function support-width
%
% Input variables
% ===============
% dim = dimension of the domain of discretization
% n = total number of nodes in the domain
% ncoord = vector of nodal coordinates 
%          1d: ncoord = [x1;x2;...;xn]
%          2d: ncoord = [x1 y1;x2 y2;...;xn yn]
%          3d: ncoord = [x1 y1 z1;x2 y2 z2;...;x3 y3 z3] 
%
% Return
% ======
% h_node = [h1; h2;...; hn]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function h_node=nodespacing(dim,n,ncoord)

clear h_node

h_node=zeros(n,1);
distances=zeros(n,n);

% STILL NAIVE IMPLEMENTATION, BUT MUCH MORE FASTER THAN IN PREVIOUS MAXENT
% VERSIONS

k=1;
for i=1:n
  % compute distances from all the nodes to node i
  distances(i,:)=sqrt(sum((ones(n,1)*ncoord(i,:)-ncoord).^2,2)).'; 
  % will avoid distance to himself being counted, otherwise distances(i,i)=0
  distances(i,i)=-1; 
  % sort the row in ascending order
  distances(i,:)=sort(distances(i,:)); % the first number in each row will be -1
  % find node spacing
  if (dim == 1)
    if (n < 3)
      h_node(i,1)=ncoord(2,1)-ncoord(1,1);
      k=k+1;
    else
      % store the first min. distances to node i (the first number is -1)
      h_node(i,1)=distances(i,2); 
      k=k+1;
    end
  elseif (dim == 2)
    % store the second min. distance to node i (the first number is -1)
    h_node(i,1)=distances(i,3); 
    k=k+1;
  elseif (dim == 3)
    % store the third min. distance to node i (the first number is -1)
    h_node(i,1)=distances(i,4); 
    k=k+1;
  else
    error('Fatal error! Dimension not yet coded.')
  end
end

clear distances



