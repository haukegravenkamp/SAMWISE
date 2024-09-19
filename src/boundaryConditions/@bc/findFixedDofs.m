function glb = findFixedDofs(obj,msh,glb)


fixedDofs = false(1,msh.nDofs);                                             % allocate fixed degrees of freedom


for i = 1:numel(obj)                                                        % loop fixed b.c.

    if isa(obj(i),'bcFixed')
        fixedDofs(obj(i).dofs) = true;                                      % set dofs to fixed
    end

end

glb.fixedDofs = fixedDofs;


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
