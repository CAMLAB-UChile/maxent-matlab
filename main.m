function main

  clc
  close all;

  % dimension of the domain of discretization
  dim = 2;         % WARNING: ilambda in prior function data must be a vector of size dim

  % number of nodes
  n = 8;

  % unstructured nodes
  unstructured=0;  % 1 = true, 0 = false
                   % WARNING: if unstructured is TRUE, then gamma in prior function
                   % data must be a vector of size number of nodes whose
                   % entries are the gamma value for each node in the mesh

  % nodal coordinates
  %ncoord = [0; 1.7857; 3.5714; 5.3571; 7.1429; 8.9286; 10.7143; 12.5000; 14.2857; 16.0714; 17.8571; 19.6429; 21.4286; 23.2143; 25.0000];
  %ncoord = [0; 0.25; 0.5; 0.75; 1];
  %ncoord = [0; 0.5; 1.0];
  %ncoord = [0; 1];
  %ncoord = [0; 0.25; 0.5; 0.75; 1];
  %ncoord = [0 0; 1 0; 1 1; 0 1];
  %ncoord = [0 0; 0.5 0; 1 0; 0 0.5; 0.5 0.5; 1 0.5; 0 1; 0.5 1; 1 1];
  %ncoord = [0.5 0.5; 1.5 0.5; 1.5 1.5; 0.5 1.5];
  %ncoord=[0 0; 1 0; 1 1; 0 1; 0.5 0.5; 0.25 0.25; 0.25 0.75; 0.75 0.25; 0.75 0.75];
  %ncoord = [0.0 0.0; 1.0 0.0; 1.0 1.0; 0.0 1.0; 0.3 0.6];
  ncoord = [0.0 0.0; 0.5 0.0; 0.0 0.5; 0.5 0.5; 1.0 0.5; 0.0 1.0; 0.5 1.0; 1.0 1.0];
  %ncoord = [0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
  %ncoord =  [0 0; 0 1.0000; 1.0000 0; 1.0000 1.0000; 0.3581  0; 0 0.5924; 0.3130 1.0000; 0.6686 0; ...
  %           1.0000 0.3339; 0.5810 0.5883; 1.0000 0.7509; 0.5983 1.0000; 0.8192 1.0000; 0.8782 0; 1.0000 0.5619; 1.0000 0.8884];

  % evaluation point
  %x = 0.297619;
  %x = 0.211325;
  %x = [0.3 0.4]; % Always Cartesian coordinates
  %x = [0.3 0.6];
  %x = [0.25 0.25];
  %x = [0.16666 0.5];
  %x=[0.5 0.5];
  x = [0.2 0.2];
  %x = [0.665465 0.447852];
  %x = [0.25 0.25];
  %x = [0.35];
  %x = [0.25 0.25 0.25];

  % prior function data
  rtol = 1e-10;                       % required tolerance for maxent iterations
  prior_type = 'gaussian';            % 'quartic_spline', 'cubic_spline', 'gaussian' or 'constant'
  
  compute = 2;                        % = 1 for basis functions computation only, 
                                      % = 2 for both basis functions and its gradient.
                                      
  gamma = 8.0*ones(size(ncoord,1),1); % support-width (dilation) parameter (needed only other than 'constant' prior)
                                      % typical values for gamma:
                                      %   gaussian = 2 to 10
                                      %   quartic spline = 1 to 2
                                      %   cubic spline = 1 to 2
                                      % Always a vector of size = number of nodes

  % initial values for lagrange multipliers
  %ilambda = 0;             % 1D
  ilambda = [0; 0];         % 2D
  %ilambda = [0; 0; 0];     % 3D

  % specify if basis functions and gradients must be printed out
  printmaxent = 'yes';  % 'yes' or 'no'

  % specify if consistency check test for phi and phider must be done
  checktest = 'yes'; % 'yes' or 'no'

  if size(x,2)~=dim
      fprintf('size(x,2) and dim dimensions mismatch: (size(x,2) = %d) vs (dim = %d)',size(x,2),dim)
  end
  if length(x)~=size(ncoord,2)
      fprintf('length(x) and size(ncoord,2) dimensions mismatch: (size(x,2) = %d) vs (size(ncoord,2) = %d)',size(x,2),size(ncoord,2))   
  end
  if size(ncoord,1)~=n
      fprintf('size(ncoord,1) and n dimensions mismatch: (size(ncoord,1) = %d) vs (n = %d)',size(ncoord,1),n)     
  end
  if length(gamma)~=n
      fprintf('length(gamma) and n dimensions mismatch: (length(gamma) = %d) vs (n = %d)',length(gamma),n)    
  end
  if length(ilambda)~=dim
      fprintf('length(ilambda) and dim dimensions mismatch: (length(ilambda) = %d) vs (dim = %d)',length(ilambda),dim)     
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %                        Compute maxent basis functions
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  % characteristic nodal spacing
  % OBS: THE RETURN OF THIS FUNCTION SHOULD BE PROVIDED FROM OUTSIDE TO MAXENT FUNCTION 
  %      SINCE IT MAY TAKE LONG TIME IF A LARGE SET OF NODES IS USED, ESPECIALLY IN 3D. BETTER IF 
  %      ANOTHER WAY TO PASS  h_node TO MAXENT FUNCTION IS EMPLOYED. FOR EXAMPLE, YOU MAY DEFINE IT 
  %      AS 2 TIMES OR 3 TIMES THE LARGEST NODAL SPACING IN YOUR MESH, AND THEN PASS IT 
  %      TO MAXENT FUNCTION AS AN ARGUMENT.
  
  h_node=nodespacing(dim,n,ncoord); % NAIVE IMPLEMENTATION, BUT MUCH MORE FASTER 
                                    % THAN PREVIOUS MAXENT VERSIONS

  maxent(dim,n,ncoord,x,prior_type,gamma,ilambda,rtol,compute,printmaxent,checktest,unstructured,h_node);

end
