function [obj,problem]=checkoptions(obj,problem)

%% CHECKOPTIONS
% Checks the validity of (some of) the options passed to samwise
% 
%
% Hauke Gravenkamp, gravenkamp.research@gmail.com 2023

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 


%% CHECK APPROACH TO SETTING MATERIAL ASSIGNMENT
checkOpt=obj.setMaterial; % current option to check
% Must be integer 0..4 or empty
if ~isempty(checkOpt)
    if ~isnumeric(checkOpt)||numel(checkOpt)~=1||~ismember(checkOpt,0:4)
        checkOpt=[];
        warning('*** Invalid value of option "setmaterial". Changing to automatic behavior.')
    end
end
if checkOpt==0&&isempty(problem.globalMaterial)
    problem.globalMaterial=1;
end
if isempty(checkOpt)
    % Automatic behavior
    if ~isempty(problem.materialFun)        % function overwrites all other definitions
        checkOpt=3;
    elseif problem.fileType==2              % material numbers from Ansys file
        checkOpt=4;
    elseif problem.fileType==3              % material numbers from image
        checkOpt=2;
    elseif ~isempty(problem.elemMaterials)  % vector of material numbers
        checkOpt=1;
    elseif ~isempty(problem.globalMaterial) % global material number
        checkOpt=0;
    else                                    % use global material number 1
        checkOpt=0;
        problem.globalMaterial=1;
        dispC('*** No valid definition of material assignment found. Choosing material 1.',obj.suppressOutput)
    end
    dispC(['*** Option "setmaterial" automatically set to ',num2str(checkOpt)],obj.suppressOutput);
    obj.setMaterial=checkOpt;
end


%% CHECK APPROACH TO SETTING ELEMENT ORDER
checkOpt=obj.setOrder;                                                      % current option to check
% Must be integer 0..3 or empty
if ~isempty(checkOpt)
    if ~isnumeric(checkOpt)||numel(checkOpt)~=1||~ismember(checkOpt,0:3)
        checkOpt=[];
        warning('*** Invalid value of option "setorder". Changing to automatic behavior.')
    end
    if checkOpt==0&&isempty(problem.globalOrder)
        problem.globalOrder=2;
        warning('*** Option "setorder" is set to 0 (global), but no global element order is defined. Choosing p=2.')
    elseif checkOpt==2&&isempty(problem.fMax)
        warning('*** Option "setorder" is set to 2 (by frequency), but no frequency fMax is defined. Choosing fMax=0.')
        problem.fMax=0;
    end
end

if isempty(checkOpt)                                                        % Automatic behavior
    if ~isempty(problem.orderFun)
        checkOpt=3;
    elseif ~isempty(problem.fMax)
        checkOpt=2;
    elseif ~isempty(problem.elemOrders)
        checkOpt=1;
    elseif ~isempty(problem.globalOrder)
        checkOpt=0;
    else
        checkOpt=0;
        problem.globalOrder=2;
        dispC('*** No valid definition of element order found. Choosing p=2.',obj.suppressOutput)
    end
    dispC(['*** Option "setorder" automatically set to ',num2str(checkOpt)],obj.suppressOutput);
    obj.setOrder=checkOpt;
end






end
