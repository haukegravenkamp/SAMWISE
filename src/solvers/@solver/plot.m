function [axisHandles, plotHandles] = plot(obj,varargin)

for i = 1:numel(obj)
    [axisHandles, plotHandles] = plotSolver(obj(i),varargin{:});
end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
