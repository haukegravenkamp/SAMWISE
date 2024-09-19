function D = computeD2d(D,assumption2D)

entriesVoigt = [1,2,6];                                                     % 2d indices in Voigt notation
T = zeros(numel(entriesVoigt),6);
for i = 1:numel(entriesVoigt)
    T(i,entriesVoigt(i)) = 1;
end
switch assumption2D

    case 'planeStrain'

        if isa(D,'function_handle')
            D = @(y,z) T*D(y,z)*T.';
        else
            D = T*D*T.';
        end

    case 'planeStress'

        Dinv = inv(D);
        D = inv(Dinv(entriesVoigt,entriesVoigt));

    case 'out-of-plane'
        D = D(entriesVoigt(end),entriesVoigt(end));

end


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
