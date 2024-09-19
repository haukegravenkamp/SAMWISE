function obj=setProperties(obj,properties)
% set properties of pde according to a list of provided options
% properties has the format {'propertyName1', value1, 'propertyName2', value2,...}

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
for i=1:floor(numel(properties)/2)                                          % loop provided options

    set(obj,properties{2*(i-1)+1},properties{2*(i-1)+2})                    % set property according to option

end
