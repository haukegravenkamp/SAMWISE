function bc = getBcLocation(obj,bc)
% get information about boundary condition

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
connectivity = obj.connectivity; 
materialNumber = obj.materialNumber;

%%
coord = obj.coord;                                                     % all vertical coordinates of vertices
allPdeDofs = [obj.pdeDofs{:}];                                              % all dofs per node for all PDEs

yMax = max(coord(:,3));
yMin = min(coord(:,3));
zMax = max(coord(:,3));
zMin = min(coord(:,3));



if strcmp(bc.location,'top')                                                % bc applied to top layer
    nodesActive = (abs(coord(:,3) - zMax)/abs(zMax)<1e-4);
elseif strcmp(bc.location,'bottom')                                         % bc applied to bottom layer
    nodesActive = (abs(coord(:,3) - zMin)/abs(zMax)<1e-4);
elseif strcmp(bc.location,'left')                                         % bc applied to bottom layer
    nodesActive = (abs(coord(:,2) - yMin)/abs(yMax)<1e-4);
elseif strcmp(bc.location,'right')                                         % bc applied to bottom layer
    nodesActive = (abs(coord(:,2) - yMax)/abs(yMax)<1e-4);
else                                                                        % undefined boundary condition
    return                                                                  % do nothing
end

nodesActive = find(nodesActive);

if isempty(bc.directions)
    bc.dofs = nonnan(allPdeDofs(nodesActive,:));                                % dofs for boundary condition (all dofs of node)
else
    bc.dofs = nonnan(allPdeDofs(nodesActive,bc.directions));                            % dofs for boundary condition (all dofs of node)
end
bc.coord = coord(nodesActive,:);                                                         % store coordinate
eleNo = any(ismember(connectivity,nodesActive),2);                         % element containing active node
bc.materialNo = materialNumber(eleNo);                                      % material number of element

end
