function [xi,eta] = triIntCoordsLocal(p,etaGLL)

nModesDom = sum(1:(p-2));                                                   % number of nodes in the interior

v = (etaGLL+1)/2;                                                                % scale to 0->1

xi  = zeros(nModesDom,1);                                                   % local coordinates on triangle (between 0 and 1!)
eta = zeros(nModesDom,1);

c = 1;                                                                      % counter
for i = 2:p                                                                 % compute local coordinates
    for j = 2:(p+1-i)
        k=p+3-i-j;
        xi(c)  = 1/3*(1 + 2*v(i) -   v(j) - v(k));
        eta(c) = 1/3*(1 -   v(i) + 2*v(j) - v(k));
        c = c+1;
    end
end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
