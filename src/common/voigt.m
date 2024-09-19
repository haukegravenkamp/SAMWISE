function x = voigt(x)

if isvector(x)

    if numel(x)==3
        x = x([1 3; 3 2]);
    else
        x = x([1 6 5; 6 2 4; 5 4 3]);
    end
else
    if numel(x)==4
        x = x([1; 4; 3]);
    else
        x = x([1, 5, 9, 8, 7, 4]);
    end

end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
