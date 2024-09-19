function [xi,yi,S]=curveNodePos(fx,df,x1,x2,etai)
%% CURVENODEPOS
% Computes the global coordinates xi,yi for GLL points along an edge.
% The edge is defined by coordinates of start and end point coordinates
% x1, x2
% Curve is defined by the function handle fx of x
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
integrand=@(x) sqrt(1 + abs(df(x)).^2);                                     % integrand (function handle)
s=@(b) integral( integrand,x1,b);                                           % arc length (integral from x1 to variable b)

S=s(x2);                                                                    % total arc length

Lx=(x2-x1);                                                                 % length in x-direction

si=(etai+1)/2*S;                                                            % arc lengths at nodal coordinates
xguess=(etai+1)/2*Lx +x1;                                                   % initial value for finding x coordinate

% find x coord of nodes corresponding to correct distribution along arc
xi=zeros(size(si));                                                         % initialize x coordinates 
for j=1:numel(si)                                                           % loop over inner nodes
    xi(j)=fzero(@(b) s(b) -  si(j), xguess(j));                             % find x coordinates
end
% yi=eval(subs(fx,xi));                                                       % corresponding y coordinates

yi = fx(xi);                                                       % corresponding y coordinates


end
