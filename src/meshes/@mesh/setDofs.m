function obj = setDofs(obj,pdes,intTab)

obj = updateConnectivity(obj,intTab);

obj = assignDofsToElements(obj,pdes);

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
