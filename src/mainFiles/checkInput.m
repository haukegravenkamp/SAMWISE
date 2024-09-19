function [geo, mat, bcd, sol, opt, res] = checkInput(inputObjects)
% check which of the essential objects required by the code are provided by
% the user. The others are created with default properties.
% If no material is provided, its name is set to 'default', meaning that
% the default material should be loaded from the database.

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
geo        = checkExistence(inputObjects, 'geometry');
bcd        = checkExistence(inputObjects, 'bc');
sol        = checkExistence(inputObjects, 'solver');
opt        = checkExistence(inputObjects, 'option');
res        = checkExistence(inputObjects, 'result');
[mat,flag] = checkExistence(inputObjects, 'material');

if ~flag
    mat.name = 'default';
    warning('no material provided, using default material parameters')
end

end


function [obj,flag] = checkExistence(inputObjects, className)

ind=find(cellfun(@(x)isa(x,className),inputObjects),1);                     % find position of given object
if isempty(ind)                                                             % if undefined
    flag = false;
    obj = eval(className);                                                  % use default
else
    flag = true;
    obj=inputObjects{ind};                                                  % use identified object
end

end
