function databasePath = getDatabasePath

%% get path to material database
currentPath = mfilename("fullpath");                                        % path to this file
materialClassPath = currentPath(1:find(currentPath=='@',1,'last')-1);       % path to parent folder of material class
databasePath = [materialClassPath,'materialDatabase'];                      % path to database
if ispc                                                                     % windows
    databasePath = [databasePath,'\'];
else                                                                        % not windows
    databasePath = [databasePath,'/']; 
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
end
