function varargout = selectMatrixIndices(dofs,varargin)

varargout = varargin;

nd = max(cellfun(@ndims,varargin));                                         % maximum number of dimensions

for i = 1:numel(varargout)

    if isempty(varargout{i})                                                % empty array in argument list
        continue                                                            % do nothing
    elseif iscell(varargout{i})                                             % argument is a cell (presumably containing arrays)
        varargout{i} = cellfun(@(x) selectMatrixIndices(dofs,x),...
            varargout{i}, 'UniformOutput', false);                          % iteratively call selectMatrixIndices for cell content
    else
        varargout{i} = varargout{i}(dofs,:);
        if nd>1
            varargout{i} = varargout{i}(:,dofs);
        end
        if nd>2
            varargout{i} = varargout{i}(:,:,dofs);
        end
        if nd>3
            varargout{i} = varargout{i}(:,:,:,dofs);
        end
    end

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
