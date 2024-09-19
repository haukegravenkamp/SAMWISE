function [coordGlobal,J,Jdet]=coordinateTransform(~,nV,coordNodes,eta12,N,Ndeta1,Ndeta2,useLinear)

if (nargin < 5)
    N = [];
    Ndeta1 = [];
    Ndeta2 = [];
    useLinear = true;
end

if nV == 3

    vert = coordNodes(1:3,2:3);
    nP = size(eta12,1);
    J = zeros(3,3,nP);
    Jdet = zeros(nP,1);
    coordGlobal = zeros(nP,3);

    for i = 1:nP
        [u,v,w,J1] = local2Barycentric(eta12(i,1),eta12(i,2));                  % compute barycentric coordinates
        [y,z,J2] = barycentric2cartesian(u,v,w,vert);
        coordGlobal(i,2:3) = [y,z];
        J(2:3,2:3,i) = J1(1:2,1:2)*J2;
        J(1,1,i) = 1;
        Jdet(i) = abs(det(J(:,:,i)));
    end

elseif nV == 4
    [coordGlobal,J,Jdet]=coordinateTransformQuad(coordNodes,eta12,N,Ndeta1,Ndeta2,useLinear);
end



end


%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
