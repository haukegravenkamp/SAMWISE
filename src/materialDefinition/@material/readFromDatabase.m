function [obj, flagRead] = readFromDatabase(obj,databasePath)
%% interface
name = obj.name;                                                            % get name of requested material

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if nargin<2
    databasePath = material.getDatabasePath;
end

%% read material from database
flagRead = false;                                                           % flag for successful reading of material
try 
    fileContent = load([databasePath,name,'.mat']);                         % try reading mat file from database
    fields = fieldnames(fileContent);                                       % variable names in file
    for i = 1:numel(fields)                                                 % loop variables
        if isa(fileContent.(fields{i}),'material')                          % if a material is found among the variables
            if flagRead                                                     % if we had found one already
                warning(['more than one material found in database...' ...  % stick with the first one, but issue warning
                    ' file, used the first one'])
                return
            end
            obj = fileContent.(fields{i});                                  % extract the variable and set it as return object
            flagRead = true;                                                % set flag
        end
    end
catch
    return
end



end
