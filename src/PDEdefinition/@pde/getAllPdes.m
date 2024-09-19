function pdes = getAllPdes(geo)
% method for initializing all PDEs
% any new PDE must be included here

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
pdes(1) = pdeAcoustics;                                                     % call acoustic PDE

if isprop(geo,'assumption2D')                                               % check if geometry includes a field for plane strain/stress assumptions
    assumption2D = geo.assumption2D;
else
    assumption2D = 'none';
end

if strcmp(assumption2D,'cylinder')
    isCylinder = true;
else
    isCylinder = false;
end



pdes(2) = pdeElasticity(assumption2D);                                      % call elastic PDE


for iPde = 1:numel(pdes)

    pdes(iPde).spatialDimension = geo.spatialDimension;
    pdes(iPde).isCylinder = isCylinder;

end



end
