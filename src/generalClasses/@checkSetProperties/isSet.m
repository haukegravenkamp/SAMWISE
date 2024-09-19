function nonEmpty = isSet(obj,properties)
% checks if properties of an object are non-empty
% properties is a cell of property names (strings)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
nProperties = numel(properties);

nonEmpty = false(nProperties,1);

for i = 1:nProperties

    nonEmpty(i) = ~isempty(obj.(properties{i}));

end

end
