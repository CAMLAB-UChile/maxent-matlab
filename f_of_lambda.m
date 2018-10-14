function [gam,dgam,hgam,phi]=f_of_lambda(dim,w,xc,lambda,len)

%
% Adapted from Marino Arroyo's matlab routines
% http://www-lacan.upc.es/arroyo/Site/Marino_Arroyo.html
%

if dim==1
  Z_I=w.*exp(-(xc*lambda));
  Z=sum(Z_I);
  phi=Z_I/Z; % contribuitng basis function evaluated at xc = xi - x
  gam=log(Z);
  temp1(1:len)=-xc.*phi;
  dgam=sum(temp1(:)); % gradient of log(Z)
  hgam = sum( phi.*(xc).^2   ) - dgam*dgam;
else
  Z_I=w.*exp(-(xc*lambda));
  Z=sum(Z_I);
  phi=Z_I/Z; % basis function
  gam=log(Z);

  for id=1:dim
    dgam(id)=sum( ( -xc(:,id) ).*phi(:) ); % gradient of log(Z)
  end
  for id=1:dim
    for jd=1:dim
        hgam(id,jd)=sum( phi(:).* ...
            ( -xc(:,id) ).*(-xc(:,jd)) )  ...
            - dgam(id)*dgam(jd);
    end
  end
end






