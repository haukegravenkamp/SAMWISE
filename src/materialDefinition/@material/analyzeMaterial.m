function obj=analyzeMaterial(obj,geo)
% check material parameters and assign missing values

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if nargin<2
    geo = [];
end

% %% get path to material database
databasePath = material.getDatabasePath;

%% set material parameters
for i=1:numel(obj)                                                          % loop over materials

    materialName = obj(i).name;                                             % read material name
    if ~isempty(materialName)                                               % material name is provided
        [obj(i), flagRead] = readFromDatabase(obj(i),databasePath);         % try to read from database
        if flagRead
            disp(['Read material "',materialName,'" from database. ' ...
                'If this is not intended, use a different name.'])
        end
    else                                                                    % no name given
        obj(i).name=['material ',num2str(i)];                               % assign name according to material number and proceed
    end


    obj(i).parameters = obj(i).(obj(i).behavior);                           % read relevant material parameters depending on type of behavior

    obj(i).parameters = assignMissingValues(obj(i).parameters,geo);         % call functions to fill parameters with missing values

    obj(i).(obj(i).behavior) = obj(i).parameters;

    if ~isempty(geo) && isprop(geo,'assumption2D') && isprop(obj(i).parameters,'D')
        obj(i).parameters.D = parametersElastic.computeD2d(obj(i).parameters.D,geo.assumption2D);
    end

end


end
