function order = getEleOrder(~,order,cMats,matNo,sizes,fMax,orderDef,orderFun)
% compute required polynomial order for an element of given length,
% frequency, and wave velocity

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
nEle = numel(sizes);

if isempty(order)
    order = zeros(nEle,1);
end

for iEle = 1:nEle                                                           % loop elements

    if order(iEle) > 0                                                      % if element order already well defined
        continue                                                            % do nothing
    end
    L = sizes(iEle);                                                        % element size
    c = min(cMats(nonnan(matNo(iEle,:))));                                  % (smallest) velocity

    if isnan(c)                                                             % wave velocity was not properly defined before
        order(iEle) = orderDef;                                             % resort to default
    else
        order(iEle) = orderFun(2*pi*fMax*L/c);                              % compute element order based on empirical function
    end

end

end
