function obj = removeNegativeAtt(obj)

kz_bottom = obj.kz_bottom;
kz_top    = obj.kz_top;
Pz_bottom = obj.Pz_bottom;
Pz_top    = obj.Pz_top;
k         = obj.k;
cp        = obj.cp;
att       = obj.att;
cg        = obj.cg;
phi       = obj.phi;
u         = obj.u;

for i = 1:size(k,2)
    indRemove = (imag(k(:,i))<1e-6);
    k(indRemove,i) = nan;                                                   % remove unwanted solutions
    cp(indRemove,i) = nan;
    att(indRemove,i) = nan;
    if ~isempty(cg)
        cg(indRemove,i) = nan;
    end
    if ~isempty(phi)
        phi(:,indRemove,i) = nan;
    end
    if ~isempty(u)
        u(:,indRemove,i) = nan;
    end
    if ~isempty(kz_bottom)
        kz_bottom(indRemove,i,:) = nan;
    end
    if ~isempty(kz_top)
        kz_top(indRemove,i,:) = nan;
    end
    if ~isempty(Pz_top)
    Pz_top(indRemove,i) = nan;
    end
    if ~isempty(Pz_bottom)
    Pz_bottom(indRemove,i) = nan;
    end
end

obj.kz_bottom = kz_bottom;
obj.kz_top    = kz_top;
obj.Pz_bottom = Pz_bottom;
obj.Pz_top    = Pz_top;
obj.k         = k;
obj.cp        = cp;
obj.att       = att;
obj.cg        = cg;
obj.phi       = phi;
obj.u         = u;

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
