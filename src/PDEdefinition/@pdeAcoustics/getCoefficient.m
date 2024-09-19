function coeffs = getCoefficient(~,mat,xyz)

c1 = mat.parameters.K;
c2 = mat.parameters.rho;

% assume constant within the element for now
c1 = c1*ones(1,1,size(xyz,1));
c2 = c2*ones(1,1,size(xyz,1));

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
coeffs = {c1,c2};

end
