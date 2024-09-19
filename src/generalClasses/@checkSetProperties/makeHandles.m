function [obj, isHandle] = makeHandles(obj,properties)
% Check if any of the properties is a function handle. If so, convert all
% of them to function handles

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

nProperties = numel(properties);

isHandle = cellfun(@(x)isa(obj.(x),'function_handle'),properties);

if any(isHandle)
    for i = 1:nProperties
        if ~isHandle(i)
            propertyValue = obj.(properties{i});
            obj.(properties{i}) = @(x,y) propertyValue;
        end
    end
    isHandle = true;

end


end
