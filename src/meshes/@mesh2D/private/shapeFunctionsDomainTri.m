function [shp,shpDeta1,shpDeta2]=shapeFunctionsDomainTri(etadi,xidi,proriolOrders,iV,indD)
%% high order shape functions on triangle

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% generate vector of basis functions and their derivatives
phi     = zeros(size(proriolOrders,1),1);                                   % basis functions
phiDxi  = zeros(size(proriolOrders,1),1);                                   % derivatives wrt xi
phiDeta = zeros(size(proriolOrders,1),1);                                   % derivatives wrt eta
for j=1:size(proriolOrders,1)
    phi(j)     = proriol(proriolOrders(j,1),proriolOrders(j,2),xidi,etadi);
    phiDxi(j)  = prorioldx1(proriolOrders(j,1),proriolOrders(j,2),xidi,etadi);
    phiDeta(j) = prorioldx2(proriolOrders(j,1),proriolOrders(j,2),xidi,etadi);
end

shp = zeros(3,size(proriolOrders,1));                                       % allocate 2D shape functions
shp(1,:)=iV*phi;                                                            % shape functions
shp(2,:)=iV*phiDxi;                                                         % derivatives wrt xi
shp(3,:)=iV*phiDeta;                                                        % derivatives wrt eta

nInt = numel(indD);
ind = size(shp,2)-nInt + (1:nInt);
shp = shp(:,ind);

shpDeta1=shp(2,:);
shpDeta2=shp(3,:);
shp=shp(1,:);

end
