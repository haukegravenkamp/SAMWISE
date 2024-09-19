function usePDE = selectPDE(~,mat)

if strcmp(mat.behavior,'acoustic')
    
    usePDE = true;
else
    usePDE = false;

end



end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
