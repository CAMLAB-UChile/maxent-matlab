function [contribute,len]=neighbors(dim,radius,xc,n)

clear contribute
node_ind=1:n; % indices of all the nodes in the domain
dist=zeros(n,1);
for i=1:dim
  dist=dist+xc(:,i).^2;
end
dist=sqrt(dist);

contribute=node_ind(dist<radius);
len=length(contribute);
clear dist

