%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            FUNCTION: maxent
%        A MATLAB implementation of the maximum-entropy basis functions
%-------------------------------------------------------------------------------
%  Version      : 3.5
%  Date         : 14-OCT-2018
%  Source code  : http://camlab.cl/software/maxent
%  Author       : A. Ortiz-Bernardin, aortizb@uchile.cl, camlab.cl/alejandro
%
%              (See Copyright and License notice in "License.txt")
%              (See updates and version details in "Version.txt")
%-------------------------------------------------------------------------------
% Purpose
% =======
% MATLAB code that computes the maximum entropy (maxent) basis functions and
% their gradients in 1D, 2D and 3D.
%
% Usage
% =====
% maxent(dim,n,ncoord,x,prior_type,gamma,ilambda,rtol,compute,printmaxent,checktest,...
%        unstructured,h_node)
%
% Input
% =====
% dim = dimension of the domain of discretization
% n = total number of nodes in the domain
% ncoord = vector of nodal coordinates
%          1d: ncoord = [x1;x2;...;xn]
%          2d: ncoord = [x1 y1;x2 y2;...;xn yn]
%          3d: ncoord = [x1 y1 z1;x2 y2 z2;...;x3 y3 z3]
% x = sampling (evaluation) point
%          1d: x = x
%          2d: x = [x y]
%          3d: x = [x y z]
% prior_type  = name of the prior function
% gamma       = parameter that controls the basis function support-width
%             = gamma in gaussian prior and
%             = alpha (dilation parameter) in cubic and quartic spline
%               priors
% ilambda = initial value of the lagrange multipliers
% checktest = "yes" to perform consistensy check on phi and phider,
%             otherwise "no"
% unstructured = "yes" or "no" (if "yes" note that it is an experimental feature. If
%                               you don't know what you are doing, I recommend using "no")
% h_node = vector containing a characteristic nodal spacing for each node in the mesh
%
% Output
% ======
%
%-------------------------------------------------------------------------------
% References
% ==========
% [1] Sukumar N. Construction of polygonal interpolants: a maximum entropy
% approach. Int. J. Numer. Meth. Engng 2004; 61(12):2159-2181.
%
% [2] Sukumar N., Wright R.W. Overview and construction of meshfree shape
% functions: From moving least squares to entropy approximants. Int. J.
% Numer. Meth. Engng 2007; 70(2):181-205.
%
% [3] Arroyo M., Ortiz M. Local maximum-entropy approximation schemes: a
% semaless bridge between finite elements and meshfree methods. Int. J.
% Numer. Meth. Engng 2006; 65(13):2167-2202.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function maxent(dim,n,ncoord,x,prior_type,gamma,ilambda,rtol,compute,printmaxent,checktest,unstructured,h_node)

% NOTE about the main variables that are computed in this function
%---------------------------------------------------------------------------------------------------
% phi = vector of maxent basis functions = [phi_1;phi_2;...;phi_len]
% phider = vector containing basis functions derivatives
%          1d: phider = [phi_1,x; phi_2,x;...; phi_len,x]
%          2d: phider = [phi_1,x phi_1,y; phi_2,x phi_2,y;...; phi_len,x phi_len,y]
%          3d: phider = [phi_1,x phi_1,y phi_1,z; phi_2,x phi_2,y phi_2,z;...; phi_len,x
%                       phi_len,y phi_len,z]
% contribute = vector containing the index of the nodes that contribute to
%              phi and phider
% len = length of contribute vector
%---------------------------------------------------------------------------------------------------

% compute basis function and its gradient
[phi,phider,contribute,len,lambda]=computephi(dim,compute,prior_type,gamma,ilambda,rtol,x,ncoord,n,h_node,unstructured);

% print phi and phider if required
if (strcmp(printmaxent,'yes'))
  if (compute == 1)
    fprintf('\n');
    fprintf('For evaluation point:\n');
    fprintf('%f\n',x);
    fprintf('Neighbor nodes: ')
    for i=1:len
      fprintf('%d ',contribute(i));
    end
    fprintf('\n');
    fprintf('MaxEnt basis functions:\n');
    for i=1:length(phi)
      fprintf('%e\n',phi(i));
    end
    fprintf('\n');   
  elseif (compute == 2)
    if (dim == 1)
      fprintf('\n')
      fprintf('For evaluation point:\n')
      fprintf('%f\n',x)
      fprintf('Neighbor nodes: ')
      for i=1:len
        fprintf('%d ',contribute(i));
      end
      fprintf('\n');
      fprintf('MaxEnt basis functions:\n')
      for i=1:length(phi)
        fprintf('%e\n',phi(i));
      end
      fprintf('\n');
      fprintf('MaxEnt basis function gradient:\n')  
      for i=1:length(phi)
        fprintf('%e\n',phider(i,1));
      end
      fprintf('\n');
    elseif (dim == 2)
      fprintf('\n');
      fprintf('For evaluation point:\n');
      fprintf('%f\n',x);
      fprintf('Neighbor nodes: ')
      for i=1:len
        fprintf('%d ',contribute(i));
      end
      fprintf('\n');
      fprintf('MaxEnt basis functions:\n');
      for i=1:length(phi)
        fprintf('%e\n',phi(i));
      end
      fprintf('MaxEnt basis function gradient:\n');
      for i=1:length(phi)
        fprintf('%e %e\n',phider(i,1),phider(i,2));
      end
      fprintf('\n');
    elseif (dim == 3)
      fprintf('\n')
      fprintf('For evaluation point:\n')
      fprintf('%f\n',x)
      fprintf('Neighbor nodes: ')
      for i=1:len
        fprintf('%d ',contribute(i));
      end
      fprintf('\n');
      fprintf('MaxEnt basis functions:\n')
      for i=1:length(phi)
        fprintf('%e\n',phi(i));
      end
      fprintf('MaxEnt basis function gradient:\n');
      for i=1:length(phi)
        fprintf('%e %e %e\n',phider(i,1),phider(i,2),phider(i,3));
      end
      fprintf('\n');
    end
  end
end

% consistency check test if specified
if (strcmp(checktest,'yes'))
  consistencycheck(dim,x,lambda,rtol,compute,phi,phider,ncoord,contribute,len)
end




