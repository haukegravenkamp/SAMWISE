function cg = computeGroupVelocity(k,u,q,M0,omega)

nK = length(k);

cg = nan(nK,1);                                                             % allocate group velocities

for j = 1:nK                                                                % loop modes
    cg(j) = -2*omega*u(:,j)'*M0*u(:,j);                                     % denominator
    cg(j) = 1i*(u(:,j)'*q(:,j)-q(:,j)'*u(:,j))/cg(j);                       % enumerator
end



end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
