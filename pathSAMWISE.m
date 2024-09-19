%% add samwise folders to path
if ispc                                                                     % windows
    addpath(genpath('.\src'))
    addpath(genpath('.\examples'))         
else                                                                        % not windows
    addpath(genpath('./src'))                
    addpath(genpath('./examples'))               
end