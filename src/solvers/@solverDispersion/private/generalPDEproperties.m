

function [cols,eleNos,dof,plotDof,variableLabel]=generalPDEproperties(pdes,cols,iPde,pde2elements,distortedGrid,acousticDisplacement,materialProvided)

if isprop(pdes(iPde),'colormap') && ~isempty(pdes(iPde).colormap)       % if PDE defines its own colormap
    cols{iPde} = pdes(iPde).colormap;                                   % use this one
end

eleNos = pde2elements{iPde};                                            % elements belonging to this PDE
dof = pdes(iPde).dof;                                                   % number of dofs per node

if isprop(pdes(iPde),'plotDof')                                         % check if dofs to be used for plotting are defined
    plotDof = pdes(iPde).plotDof;                                       % use defined dofs
else
    plotDof = 1:dof;                                                    % use all dofs given by pde
end

if (distortedGrid || acousticDisplacement) ...
        &&  isa(pdes(iPde),'pdeAcoustics')                               % show distorted mesh or plot displacements in acoustics problem
    if materialProvided                                                 % check if material is proviopt.plotting.ded (required for displacement computation)
        computeDisplacement = true;                                     % flag for computing displacements
    else
        warning(['To compute displacements for acoustic wave ' ...
            'equation material objects must be passed to ' ...
            'plotWavefield. Option distortedGrid will be ignored ' ...
            'for acoustic PDE.'])
        computeDisplacement = false;                                    % don't compute displacement if material is missing
    end
else
    computeDisplacement = false;
end

if acousticDisplacement && computeDisplacement                          % displacements to be plotted in acoustic problem
    variableLabel = '$u$';                                              % always call this unknown 'u'
else                                                                    % all other cases
    variableLabel = pdes(iPde).variableName;                            % use variable name as defined by PDE
end


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
