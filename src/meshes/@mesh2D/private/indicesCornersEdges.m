function indE = indicesCornersEdges(p,nC)

%% get indices of corner and edge modes
startEdge = cumsum([0,p(3:end)-1])+nC;
indE = nan(nC,max(p(3:end))+1);
for iSide = 1:nC
    ind = [iSide,startEdge(iSide)+(1:(p(2+iSide)-1)), 1 + mod(iSide, nC)];
    indE(iSide,1:numel(ind)) = ind;
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
end
