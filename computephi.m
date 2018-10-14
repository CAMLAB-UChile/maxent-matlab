function [phi,phider,contribute,len,lambda]=computephi(dim,compute,prior_type,gamma,ilambda,rtol,x,ncoord,n,h_node,unstructured)

clear phi
clear phider
clear contribute
clear len

% compute prior function 
[w,wder,xc_c,contribute,len]=prior(dim,prior_type,gamma,rtol,x,ncoord,n,h_node,unstructured);

lambda=ilambda;
R=10*ones(1,dim);
niter=0;

%Newton iteration: Adapted from Marino Arroyo's matlab routines
%(http://www-lacan.upc.es/arroyo/Site/Marino_Arroyo.html)
dlam=10*ones(1,dim);
while (norm(R)>rtol)
  [gam,R,J,phi]=f_of_lambda(dim,w,xc_c,lambda,len); % f(lambda) = Grad(log Z) = 0
  if (abs(rcond(J)))<1e-8
    disp('Newton Failed, near to singular matrix')
  end
  dlam=-J\R';
  lambda=lambda+dlam;
  niter=niter+1;
  if (niter>100) 
    disp('Newton Failed 2, no convergence in 100 iterations')
  end
end
lambda=lambda-dlam;

fprintf('Newton interations: %d \n',niter)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute phider in an matrix nx2 as follows: phider = [ ... ... ; der_x  der_y; ... ...]
%
% General form of phider (example: 2D)
%
%      = [phi phi].*(xc_c*(invH-invH*A))+[phi phi].*MA-phi*MC
%
% this is arranged to compute phider in the format [ ... ... ; der_x  der_y; ... ...]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (compute == 1)
  phider=[];
elseif (compute == 2)
  % compute hessian and matrices to compute derivatives for any prior
  [H,MA,A,MC]=phidermatrices(dim,w,wder,len,xc_c,phi); 
  % compute basis functions derivatives
  if (dim == 1) 
    invH=1/H(1,1);
    phider=phi.*(xc_c*(invH-invH*A))+phi.*MA-phi*MC;
  elseif (dim == 2)      
    invH=inv(H); 
    phider=[phi phi].*(xc_c*(invH-invH*A))+[phi phi].*MA-phi*MC;
  elseif (dim == 3)
    invH=inv(H);
    phider=[phi phi phi].*(xc_c*(invH-invH*A))+[phi phi phi].*MA-phi*MC;
  else
    error('Fatal error! Dimension not yet coded.')
  end
end
