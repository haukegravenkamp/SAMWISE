function [sol, bcd, glb] = callSolver(obj,glb,msh,geo,bcd,mat,opt)

multiWaitbar('CloseAll');                                                   % close previous waitbars if any

sol(numel(glb)) = solver;
nSet = numel(glb);

for i = 1:nSet                                                              % loop global systems, e.g. different circumferential orders
    [sol(i), bcd, glb(i)] = computeSolution(obj,glb(i),msh,geo,bcd,mat,opt);
end

multiWaitbar('CloseAll');                                                   % close all waitbars

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
