function [obj,fMax] = assignElementOrder(obj,geo,mat,opt,fMax)
% assign interpolation order to elements
% affects only elements that do not have a valid element order already
% defined

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
nEle = obj.nEle;                                                            % number of elements

eleOrderDefault = opt.numerics.eleOrderDefault;
eleOrderFunction = opt.numerics.eleOrderFunction;

eleSizes = obj.size;
materialNumber = obj.materialNumber;

usedMaterials = [geo.layers.material];

eleOrder = obj.eleOrder;
if numel(eleOrder)<nEle
    eleOrder = [eleOrder;zeros(nEle-numel(eleOrder),1)];
end

%% get wave velocity relevant for computing dimensionless frequency
cMats = getWaveVelocity(mat,usedMaterials);
fMax = obj.getFmax(fMax,eleSizes,cMats,materialNumber);
eleOrder = getEleOrder(obj,eleOrder,cMats,materialNumber,eleSizes,fMax,eleOrderDefault,eleOrderFunction);

%% output
obj.eleOrder = eleOrder;
obj.eleOrderAll = eleOrder;
obj.nEleNodes = eleOrder+1;

end
