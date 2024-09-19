function [shp,shpDeta1,shpDeta2]=interiorModesTri(shp,shpDeta1,shpDeta2,eta1,eta2,indD,p,proriolOrders,iV)


%% add interior modes

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if p>2

    % The inner modes are based on the standard triangular elements and use
    % a different coordinate transformation. Hence the following overly
    % complicated back and forth transformations

    [u,v,w,J1] = local2Barycentric(eta1,eta2);                              % compute barycentric coordinates
    vert = [0 0; 1 0; 0 1];
    [xi,etaTri,J2] = barycentric2cartesian(u,v,w,vert);                     % used here to map barycentric to local triangle coordinates
    [xid,etad,J3]  = tri2square(xi,etaTri);                                 % map on square

    [shpDom,shpDeta1Dom,shpDeta2Dom]=shapeFunctionsDomainTri(etad,xid,proriolOrders,iV,indD); % compute shape functions for inner nodes

    shp(indD)=shpDom;
    shpDxy=J1(1:2,1:2)*J2*J3*[shpDeta1Dom;shpDeta2Dom];
    shpDeta1(indD)=shpDxy(1,:);
    shpDeta2(indD)=shpDxy(2,:);

end

end
