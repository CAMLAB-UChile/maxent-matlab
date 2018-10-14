function consistencycheck(dim,x,lambda,rtol,compute,phi,phider,coord,contribute,len)

num=round(-log10(rtol));
num=round(num/2.0);
tol=10^-(num);

if (dim == 1)
  xcoord=zeros(len,1); % store the x-coordinate of the contributing nodes
  for i=1:len
    xcoord(i,1)=coord(contribute(i),1);
  end

  if (compute == 1)
    sumphi=sum(phi);
    sumphixi=dot(phi,xcoord);
    if ((abs(sumphi-1) < tol) && (abs(sumphixi-x(1)) < tol))
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = %f\n',x)
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('\n')
    end  
  elseif (compute == 2)
    phiderx(1:len,1)=phider(1:len,1); % vector containing dphi/dx for each of the contributing nodes
    sumphi=sum(phi);
    sumphixi=dot(phi,xcoord);
    sumphiderxxi=dot(phiderx,xcoord);
    if ( (abs(sumphi-1) < tol) && (abs(sumphixi-x(1)) < tol) &&...
         (abs(sumphiderxxi-1) < tol) )
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = %f\n',x)
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('sum(phiderx*xi) = %f\n',sumphiderxxi)
      fprintf('\n')
    end   
  end
elseif (dim == 2)
  xcoord=zeros(len,1); % store the x-coordinate of the contributing nodes
  ycoord=zeros(len,1); % store the y-coordinate of the contributing nodes
  for i=1:len
    xcoord(i,1)=coord(contribute(i),1);
    ycoord(i,1)=coord(contribute(i),2);
  end

  if (compute == 1)
    sumphi=sum(phi);
    sumphixi=dot(phi,xcoord);
    sumphiyi=dot(phi,ycoord);
    if ((abs(sumphi-1.0) < tol) && (abs(sumphixi-x(1)) < tol) &&...
        (abs(sumphiyi-x(2)) < tol))
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = (%f,%f)\n',x(1),x(2))
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('sum(phi*yi) = %f\n',sumphiyi)   
      fprintf('\n')
    end   
  elseif (compute == 2)
    phiderx(1:len,1)=phider(1:len,1); % vector containing dphi/dx for each of the contributing nodes
    phidery(1:len,1)=phider(1:len,2); % vector containing dphi/dy for each of the contributing nodes
    sumphi=sum(phi);
    sumphixi=dot(phi,xcoord);
    sumphiyi=dot(phi,ycoord);
    sumphiderxxi=dot(phiderx,xcoord);
    sumphideryxi=dot(phidery,xcoord);
    sumphiderxyi=dot(phiderx,ycoord);
    sumphideryyi=dot(phidery,ycoord);
    if ((abs(sumphi-1.0) < tol) && (abs(sumphixi-x(1)) < tol) &&...
        (abs(sumphiyi-x(2)) < tol) &&  (abs(sumphiderxxi-1.0) < tol) &&...
        (abs(sumphideryxi-0.0) < tol) && (abs(sumphiderxyi-0.0) < tol) &&...
        (abs(sumphideryyi-1.0) < tol))
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = (%f,%f)\n',x(1),x(2))
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('sum(phi*yi) = %f\n',sumphiyi)   
      fprintf('sum(phiderx*xi)  sum(phidery*xi) = %f  %f\n',sumphiderxxi,sumphideryxi)
      fprintf('sum(phiderx*yi)  sum(phidery*yi) = %f  %f\n',sumphiderxyi,sumphideryyi)     
      fprintf('\n')
    end    
  end
elseif (dim == 3)
  xcoord=zeros(len,1); % store the x-coordinate of the contributing nodes
  ycoord=zeros(len,1); % store the y-coordinate of the contributing nodes
  zcoord=zeros(len,1); % store the z-coordinate of the contributing nodes
  for i=1:len
    xcoord(i,1)=coord(contribute(i),1);
    ycoord(i,1)=coord(contribute(i),2);
    zcoord(i,1)=coord(contribute(i),3);
  end
  
  if (compute == 1)
    sumphi = sum(phi);
    sumphixi=dot(phi,xcoord);
    sumphiyi=dot(phi,ycoord);
    sumphizi=dot(phi,zcoord);
    if ((abs(sumphi-1.0) < tol) && (abs(sumphixi-x(1)) < tol) &&...
        (abs(sumphiyi-x(2)) < tol) && (abs(sumphizi-x(3)) < tol))
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = (%f,%f,%f)\n',x(1),x(2),x(3))
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('sum(phi*yi) = %f\n',sumphiyi)
      fprintf('sum(phi*zi) = %f\n',sumphizi)
      fprintf('\n')
    end     
  elseif (compute == 2)
    phiderx(1:len,1)=phider(1:len,1); % vector containing dphi/dx for each of the contributing nodes
    phidery(1:len,1)=phider(1:len,2); % vector containing dphi/dy for each of the contributing nodes
    phiderz(1:len,1)=phider(1:len,3); % vector containing dphi/dz for each of the contributing nodes
    sumphi = sum(phi);
    sumphixi=dot(phi,xcoord);
    sumphiyi=dot(phi,ycoord);
    sumphizi=dot(phi,zcoord);
    sumphiderxxi=dot(phiderx,xcoord);
    sumphideryxi=dot(phidery,xcoord);
    sumphiderzxi=dot(phiderz,xcoord);
    sumphiderxyi=dot(phiderx,ycoord);
    sumphideryyi=dot(phidery,ycoord);
    sumphiderzyi=dot(phiderz,ycoord);
    sumphiderxzi=dot(phiderx,zcoord);
    sumphideryzi=dot(phidery,zcoord);
    sumphiderzzi=dot(phiderz,zcoord);
    if ((abs(sumphi-1.0) < tol) && (abs(sumphixi-x(1)) < tol) &&...
        (abs(sumphiyi-x(2)) < tol) && (abs(sumphizi-x(3)) < tol) &&...
        (abs(sumphiderxxi-1.0) < tol) &&...
        (abs(sumphideryxi-0.0) < tol) &&...
        (abs(sumphiderzxi-0.0) < tol) &&...
        (abs(sumphiderxyi-0.0) < tol) &&...
        (abs(sumphideryyi-1.0) < tol) &&...
        (abs(sumphiderzyi-0.0) < tol) &&...
        (abs(sumphiderxzi-0.0) < tol) &&...
        (abs(sumphideryzi-0.0) < tol) &&...
        (abs(sumphiderzzi-1.0) < tol))
      fprintf('Consistency check ... OK\n')
    else
      fprintf('\n')
      fprintf('For evaluation point = (%f,%f,%f)\n',x(1),x(2),x(3))
      fprintf('Consistency check ... FAILED\n')
      fprintf('Lagrange multiplier = %f\n',lambda)
      fprintf('sum(phi) = %f\n',sumphi)
      fprintf('sum(phi*xi) = %f\n',sumphixi)
      fprintf('sum(phi*yi) = %f\n',sumphiyi)
      fprintf('sum(phi*zi) = %f\n',sumphizi)
      fprintf('sum(phiderx*xi)  sum(phidery*xi)  sum(phiderz*xi) = %f  %f  %f\n',sumphiderxxi,sumphideryxi,sumphiderzxi)
      fprintf('sum(phiderx*yi)  sum(phidery*yi)  sum(phiderz*yi) = %f  %f  %f\n',sumphiderxyi,sumphideryyi,sumphiderzyi)
      fprintf('sum(phiderx*zi)  sum(phidery*zi)  sum(phiderz*zi) = %f  %f  %f\n',sumphiderxzi,sumphideryzi,sumphiderzzi)
      fprintf('\n')
    end 
  end
end
