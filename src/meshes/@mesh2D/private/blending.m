function iT = blending(p,coord,nModes,nDom,intTab)

T=eye(nModes);

modesStart = nModes-nDom;
for i=1:nDom                                                                % loop inner modes

    shp = shapeFunctionsXnY(p,coord(i,:),nModes,nDom,intTab,[]);
    T(modesStart+i,:) = shp(1,:);

end

iT=inv(T);

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
