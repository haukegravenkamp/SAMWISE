function geo = hRefine(geo,opt)

hRef = opt.numerics.hRefinement;                                            % requested number of refinement steps

if hRef == 0                                                                % no refinement
    return                                                                  % do nothing
end

for i=1:numel(geo.layers)                                                   % loop layers

    geo.layers(i).nEle = geo.layers(i).nEle*2^hRef;                         % refine

end



end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
