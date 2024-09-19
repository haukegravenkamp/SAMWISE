function coeffs = getCoefficient(~,mat,xyz)

c1 = mat.parameters.D;
c2 = mat.parameters.rho;

if isa(c1,'function_handle')
    c1Ele = zeros(size(c1(0,0),1),size(c1(0,0),2),size(xyz,1));
    for i = 1:size(xyz,1)
        c1Ele(:,:,i) = c1(xyz(i,2),xyz(i,3));
    end
else                                                                        % assume constant within the element
    c1Ele = repmat(c1,1,1,size(xyz,1));
end

if isa(c2,'function_handle')
    c2Ele = zeros(size(c2,1),size(c2,2),size(xyz,1));
    for i = 1:size(xyz,1)
        c2Ele(:,:,i) = c2(xyz(i,2),xyz(i,3));
    end
else                                                                        % assume constant within the element
    c2Ele = c2*ones(1,1,size(xyz,1));
end
coeffs = {c1Ele,c2Ele};

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
