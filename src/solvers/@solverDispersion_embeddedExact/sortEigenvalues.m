function obj = sortEigenvalues(obj)

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
    [~,ind] = sort(real(k(:,i)));
    k(:,i) = k(ind,i);                                                      % remove unwanted solutions
    cp(:,i) = cp(ind,i);
    att(:,i) = att(ind,i);
    if ~isempty(cg)
        cg(:,i) = cg(ind,i);
    end
    if ~isempty(phi)
        phi(:,:,i) = phi(:,ind,i);
    end
    if ~isempty(u)
        u(:,:,i) = u(:,ind,i);
    end
    if ~isempty(kz_bottom)
        kz_bottom(:,i,:) = kz_bottom(ind,i,:);
    end
    if ~isempty(kz_top)
        kz_top(:,i,:) = kz_top(ind,i,:);
    end
    if ~isempty(Pz_top)
        Pz_top(:,i) = Pz_top(ind,i);
    end
    if ~isempty(Pz_bottom)
        Pz_bottom(:,i) = Pz_bottom(ind,i);
    end

end

%% remove modes that only have nans
indExist = any(~isnan(k),2);
k = k(indExist,:);
att = att(indExist,:);
cp = cp(indExist,:);
if ~isempty(cg)
    cg = cg(indExist,:);
end
if ~isempty(phi)
    phi = phi(:,indExist,:);
end
if ~isempty(u)
    u = u(:,indExist,:);
end
if ~isempty(kz_bottom)
    kz_bottom = kz_bottom(indExist,:,:);
end
if ~isempty(kz_top)
    kz_top = kz_top(indExist,:,:);
end
if ~isempty(Pz_top)
    Pz_top = Pz_top(indExist,:);
end
if ~isempty(Pz_bottom)
    Pz_bottom = Pz_bottom(indExist,:);
end

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% output
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
