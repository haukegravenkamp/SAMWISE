function color=myColor(i,bright)
% A few selected colors that I like to use for plots

if nargin<2
    bright = true;
end



% brighter variant to better distinguish in paper
if bright
colors=[...
    255*[0.44 0.44 0.44];...               % black
    255*[0.89 0.26 0.30];...            % red
    43 125 161;...        % blue
   255*[0.05 0.62 0.31];...             % green
    150 30 180;...          % purple
    215 85 25;...           % orange
    128 128 128]/255;       % gray

else
colors=[...
    0 0 0;...               % black
    162 20 47;...           % red
    % 0 66 160;...            % blue
    43, 125, 161;...
    0 80 36;...             % green
    150 30 180;...          % purple
    215 85 25;...           % orange
    128 128 128]/255;       % gray

end



ind=mod(i-1,size(colors,1))+1;
color=colors(ind,:);
