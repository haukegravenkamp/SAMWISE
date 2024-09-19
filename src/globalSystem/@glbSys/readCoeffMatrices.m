function coeffM = readCoeffMatrices(obj)

E0  = obj.feMatrices{1,1};
E11 = obj.feMatrices{2,1};
E12 = obj.feMatrices{3,1};
E2  = obj.feMatrices{4,1};
C0  = obj.feMatrices{1,2};
if isempty(C0)
    C0 = 0*E0;
end
M0  = obj.feMatrices{1,3};
if isempty(M0)
    M0 = 0*E0;
end

coeffM = {E0,E11,E12,E2,C0,M0};

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
