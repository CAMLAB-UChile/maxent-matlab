% function [w,wder,xc,contribute,len] = 
%                prior(dim, prior_type,gamma,h_node,xc,n)
%
% Purpose
% =======
% Compute prior functions (wc) and their derivatives (wderc), remove nodal 
% information of those nodes that do not contribute to phi(x), construct
% contribute vector that contains the index of nodes contributing to
% phi(x) and its lenght (r)
%
% Input variables
% ===============
% dim         = dimension of the domain of discretization
% prior_type  = name of the prior function
% gamma       = parameter that controls the basis function support-width
%             = gamma in gaussian prior and 
%             = alpha (dilation parameter) in cubic and quartic spline
%               priors
%             = if unstructured = false, then gamma is a unique value for
%               all the nodes in the domain. Otherwise, it must be a vector
%               containing a value for each node
% rtol        = required newton tolerance (used to compute Gaussian prior
%               support radius)
% h_node      = vector containing characteristic nodal space
%             = [h1; h2;...; hn]                                
% x           = sampling (evaluation) point
%               1d: x = x
%               2d: x = [x y]
%               3d: x = [x y z]
% n           = total number of nodes in the domain
%
% Return
% ======
% wc         = [w1; w2;...; wn]
% 1d: wcder  = [w1,x; w2,x;...;wn,x]
% 2d: wcder  = [w1,x w1,y; w2,x w2,y;...;wn,x wn,y]
% 3d: wcder  = [w1,x w1,y w1,z; w2,x w2,y w2,z;...;wn,x wn,y wn,z]
% xc         = shifted coordinates containing the information only of the
%              nodes that contributes to phi(x)
% contribute = [node_1;...;node_m], m = number of contributing nodes
% len        = length of contribute vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [w,wder,xc,contribute,len] = prior(dim,prior_type,gamma,rtol,x,ncoord,n,h_node,unstructured)

clear contribute
clear w
clear wder
clear xc

xc=shiftcoordinates(dim,n,ncoord,x); % shifted coordinates

if ( strcmp(prior_type,'constant') )
  w=ones(n,1);
  wder=zeros(n,dim);
  contribute=1:n;
  len=n;
elseif ( strcmp(prior_type,'gaussian') )
  
  if unstructured % each node with their own beta
    tolzero=rtol; % tolerance at which Gaussian prior is considered zero
    radius=zeros(n,1);
    for i=1:n
      radius(i)=h_node(i)*sqrt(-log(tolzero)/gamma(i)); 
    end
    [contribute,len]=neighbors(dim,radius,xc,n);
    beta=zeros(len,1);
    for j=1:len
      beta(j)=gamma(contribute(j))/(h_node(contribute(j))^2);
    end  
  else % h constant but each node with their own gamma and therefore their own beta
    h=max(h_node);
    tolzero=rtol; % tolerance at which Gaussian prior is considered zero
    radius=zeros(n,1);
    for i=1:n
      radius(i)=h*sqrt(-log(tolzero)/gamma(i)); 
    end
    [contribute,len]=neighbors(dim,radius,xc,n);
    beta=zeros(len,1);
    for j=1:len
      beta(j)=gamma(contribute(j))/(h*h);
    end  
  end
    
  % compute Gaussian prior and its gradient
  if (dim == 1)
    xc=xc(contribute);
    w=exp(-beta.*(xc.^2));
    wder=2*(w.*xc).*beta;    
  elseif (dim == 2)
    xc=xc(contribute,:);
    w=exp(-beta.*(sum((xc.').^2).')); 
    wder=2*([w w].*xc).*[beta beta]; 
  elseif (dim == 3)
    xc=xc(contribute,:);
    w=exp(-beta.*(sum((xc.').^2).')); 
    wder=2*([w w w].*xc).*[beta beta beta]; 
  else
    error('Fatal error! Dimension not yet coded.')
  end
elseif ( strcmp(prior_type,'quartic_spline') )
  
  if unstructured % each node with their own h and their own beta and therefore with their own dmi = 1/radius_of_basis_function
    radius=gamma.*h_node; 
    [contribute,len]=neighbors(dim,radius,xc,n);
    dmi=zeros(len,1);
    for j=1:len
      dmi(j)=1/(gamma(contribute(j))*h_node(contribute(j)));
    end 
  else % h constant but each node with their own gamma and therefore with their own dmi = 1/radius_of_basis_function
    h=max(h_node);
    radius=gamma.*h;
    [contribute,len]=neighbors(dim,radius,xc,n);
    dmi=zeros(len,1);
    for j=1:len
      dmi(j)=1/(gamma(contribute(j))*h);
    end 
  end
  
  % compute main variable of the quartic_spline
  if (dim == 1)
    xc=xc(contribute);
    ri=dmi.*sqrt(xc.^2);       
  else
    xc=xc(contribute,:);
    ri=dmi.*sqrt(sum((xc.').^2).');
  end

  % compute quartic spline prior and its gradient
  w=1-6*ri.^2+8*ri.^3-3*ri.^4; 
  const=12*(dmi.^2)-(24*ri).*(dmi.^2)+(12*(ri.^2)).*(dmi.^2);
  if (dim == 1)
    wder=const.*xc;
  elseif (dim == 2)
    wder=[const const].*xc; 
  elseif (dim == 3)
    wder=[const const const].*xc; 
  end    
elseif ( strcmp(prior_type,'cubic_spline') )

  if unstructured % each node with their own h and their own beta and therefore with their own dmi = 1/radius_of_basis_function
    radius=gamma.*h_node; 
    [contribute,len]=neighbors(dim,radius,xc,n);
    dmi=zeros(len,1);
    for j=1:len
      dmi(j)=1/(gamma(contribute(j))*h_node(contribute(j)));
    end 
  else % h constant but each node with their own gamma and therefore with their own dmi = 1/radius_of_basis_function
    h=max(h_node);
    radius=gamma.*h;
    [contribute,len]=neighbors(dim,radius,xc,n);
    dmi=zeros(len,1);
    for j=1:len
      dmi(j)=1/(gamma(contribute(j))*h);
    end 
  end  

  % compute main variable of the cubic_spline
  if (dim == 1)
    xc=xc(contribute);
    ri=dmi.*sqrt(xc.^2);       
  else
    xc=xc(contribute,:);
    ri=dmi.*sqrt(sum((xc.').^2).');
  end

  % compute cubic spline prior and its gradient
  w=zeros(len,1);
  wder=zeros(len,1);
  for i = 1:len
    if (ri(i,1) <= 0.5)
      w(i,1)=2/3-4*ri(i,1)^2+4*ri(i,1)^3; % w = [ w1 ; w2 ; ... ; wn]
      const=8-12*ri(i,1);
      if (dim == 1)
        wder(i,1)=( (dmi(i,1)^2)*const ).*xc(i,1);
      elseif (dim == 2)
        wder(i,1:2)=( (dmi(i,1)^2)*[const const] ).*xc(i,1:2);
      elseif (dim == 3)
        wder(i,1:3)=( (dmi(i,1)^2)*[const const const] ).*xc(i,1:3);
      end
    else
      w(i,1)=4/3-4*ri(i,1)+4*ri(i,1)^2-4*(ri(i,1)^3)/3; % w = [ w1 ; w2 ; ... ; wn]       
      const=4/ri(i,1)-8+4*ri(i,1);
      if (dim == 1)
        wder(i,1)=( (dmi(i,1)^2)*const ).*xc(i,1);
      elseif (dim == 2)
        wder(i,1:2)=( (dmi(i,1)^2)*[const const] ).*xc(i,1:2);
      elseif (dim == 3)
        wder(i,1:3)=( (dmi(i,1)^2)*[const const const] ).*xc(i,1:3);
      end
    end
  end
end

% check if there are enough neighbors to compute MaxEnt basis functions

%{
if ( (dim == 1) && (len < 2) )
  if (strcmp(prior_type,'gaussian'))
      error('MaxEnt solution will not converge ... need more neighbors. Solution: decrease parameter gamma.')
  else
      error('MaxEnt solution will not converge ... need more neighbors. Solution: increase parameter gamma.')
  end
elseif ( (dim == 2) && (len < 3) )
  if (strcmp(prior_type,'gaussian'))
      error('MaxEnt solution will not converge ... need more neighbors. Solution: decrease parameter gamma.')
  else
      error('MaxEnt solution will not converge ... need more neighbors. Solution: increase parameter gamma.')
  end
elseif ( (dim == 3) && (len < 4) )
  if (strcmp(prior_type,'gaussian'))
      error('MaxEnt solution will not converge ... need more neighbors. Solution: decrease parameter gamma.')
  else
      error('MaxEnt solution will not converge ... need more neighbors. Solution: increase parameter gamma.')
  end
end
%}


  


