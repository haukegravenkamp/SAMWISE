function col=myBlues(n,perc)
% colors based on some palette found online
% b are some base color
% if n is provided, base colors are interpolated to n colors
% perc returns only frequencies between perc(1) and perc(2) percent of the
% full range

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
b=[ ...
1, 42, 74;...
1, 58, 99;...
1, 73, 124;...
1, 79, 134;...
42, 111, 151;...
44, 125, 160;...
70, 143, 175;...
97, 165, 194;...
137, 194, 217;...
169, 214, 229 ...
]/255;                                                                      % base colors

if (nargin <1) || isempty(n)                                                % no interpolation
    col=b;                                                                  % return base colors
else                                                                        % interpolate to n points
    nB=size(b,1);
    b1=interp1(linspace(1,n,nB)',b(:,1),(1:n)');
    b2=interp1(linspace(1,n,nB)',b(:,2),(1:n)');
    b3=interp1(linspace(1,n,nB)',b(:,3),(1:n)');
    col=[b1,b2,b3];
end

if (nargin > 1)  && ~isempty(perc)                                          % remove percentages
    nC=size(col,1);
    ind = round(perc(1)*nC/100):round(perc(2)*nC/100);
    col=col(ind,:);
    

end

end
