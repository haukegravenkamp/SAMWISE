function obj = getDerivatives(obj,msh,glb)

obj.uDy{msh.nEle} = {};                                                     % allocate derivatives for each elements

u = obj.u;                                                                  % solution
dofPlate = unique(nonnan(msh.dofsGlobal(:)));                               % all dofs of plate

for iEle = 1:msh.nEle                                                       % loop elements

    shpDy = glb.shpDy{iEle};                                                % read shape function derivatives
    if isempty(shpDy)
        continue
    end
    uDyC = zeros(size(u));                                                  % allocate derivatives of current element

    for i = 1:size(uDyC,2)                                                  % loop modes

        for j = 1:size(uDyC,3)                                              % loop frequencies
            uDyC(dofPlate,i,j) = shpDy*u(dofPlate,i,j);
        end
    end
    obj.uDy{iEle} = uDyC;
end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
