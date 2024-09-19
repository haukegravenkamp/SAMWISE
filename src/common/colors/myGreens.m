function col=myGreens(n,perc)
% just a modification of the blue colormaps

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if (nargin < 1) 
    n = [];
end
if (nargin < 2) 
    perc = [];
end

col=myBlues(n,perc);
col = col(:,[1,3,2])*0.85;

end
