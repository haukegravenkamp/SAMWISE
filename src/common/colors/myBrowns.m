function col=myBrowns(n,perc)
% colors based on some palette found online
% https://coolors.co/palette/fcdeba-e7c7a5-d2b18f-bd9a7a-a98365-946c4f-7f563a-6a3f24-55280f
% b are some base color
% if n is provided, base colors are interpolated to n colors
% perc returns only frequencies between perc(1) and perc(2) percent of the
% full range

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
b=[ ...
85, 40, 15;...
106, 63, 36;...
127, 86, 58;...
148, 108, 79;...
169, 131, 101;...
189, 154, 122;...
210, 177, 143;...
231, 199, 165;...
252, 222, 186;...
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
