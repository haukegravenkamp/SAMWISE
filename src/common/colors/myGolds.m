function col=myGolds(~,perc)
% this is just the gray colorbar multiplied by a brownish color

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
col = (gray.*[189, 154, 122]+30)/255;

if (nargin > 1)  && ~isempty(perc)                                          % remove percentages
    nC=size(col,1);
    ind = round(perc(1)*nC/100):round(perc(2)*nC/100);
    col=col(ind,:);
end

end
