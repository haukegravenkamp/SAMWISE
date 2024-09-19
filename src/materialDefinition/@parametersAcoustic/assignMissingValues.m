
function obj = assignMissingValues(obj,~)
% standard acoustic material: bulk modulus K, wave velocity cl, density rho
% two of the constants must be given, third one is computed

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if isSet(obj,{'K','rho'})                                                   % given K, rho
    obj.cl = sqrt(obj.K/obj.rho);                                           % compute cl
elseif isSet(obj,{'cl','rho'})                                              % given cl, rho
    obj.K = obj.cl^2*obj.rho;                                               % compute K
elseif isSet(obj,{'K','cl'})                                                % given K, cl
    obj.rho = obj.K/obj.cl^2;                                               % compute rho
else
    error('two out of three constants required for acoustic material')
end
end
