function OneDshp = shpLagrange(nNodes,x,IntgTable)
%% SHPLAGRANGE
% Lagrange shape functions
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

OneDshp = zeros(2,nNodes);
OneDshp(1,:) = 1;

pts = IntgTable(nNodes).GLL.xi;
prd = IntgTable(nNodes).GLL.prd;

%% efficient vectorized version
a=~eye(nNodes).*(x-pts')+eye(nNodes);
OneDshp(1,:)=prod(a).*prd;
ind=(a==0);
nz=all(~ind);
OneDshp(2,nz) =OneDshp(1,nz).*(sum(1./a(:,nz))-1);
a(ind)=1;
OneDshp(2,~nz)=prod(a(:,~nz)).*prd(~nz);


end

