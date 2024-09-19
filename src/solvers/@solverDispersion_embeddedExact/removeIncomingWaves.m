function obj = removeIncomingWaves(obj,opt)

%% interface
kz_bot    = obj.kz_bottom;
kz_top    = obj.kz_top;
Pz_bot    = obj.Pz_bottom;
Pz_top    = obj.Pz_top;
k         = obj.k;
cp        = obj.cp;
att       = obj.att;
cg        = obj.cg;
phi       = obj.phi;
u         = obj.u;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
unboundedUsePoyntingVec = opt.numerics.unboundedUsePoyntingVec;
unboundedRemoveThreshold = opt.numerics.unboundedRemoveThreshold;

%% check which variables exist
if isempty(Pz_bot)
    existPB = false;
else
    existPB = true;
end
if isempty(Pz_top)
    existPT = false;
else
    existPT = true;
end
if isempty(kz_bot)
    existkB = false;
else
    existkB = true;
end
if isempty(kz_top)
    existkT = false;
else
    existkT = true;
end
if isempty(cg)
    existCg = false;
else
    existCg = true;
end
if isempty(phi)
    existPhi = false;
else
    existPhi = true;
end
if isempty(u)
    existU = false;
else
    existU = true;
end


maxKunb = max([...
    max(abs( nonnan(kz_bot(:)) )) , ...
    max(abs( nonnan(kz_top(:)) )) ]);                                       % maximum of all wavenumbers in unbounded domains
thresK = unboundedRemoveThreshold*maxKunb;                                  % threshold for determining direction of wavenumber

% maxP = max(max(abs(Py_bot(:))),max(abs(Py_top(:))));
maxP = max(abs([Pz_bot(:);Pz_top(:)]));
thresP = unboundedRemoveThreshold*maxP;                                  % threshold for determining direction of wavenumber


sizeK = size(k,1);
indB = true(sizeK,1);
indT = true(sizeK,1);

for i = 1:size(k,2)                                                         % loop frequencies
    if unboundedUsePoyntingVec
        if existPB                                                          % power flux at bottom surface defined
            indB = Pz_bot(:,i)<-thresP;                                     % negative power flux at bottom (away from plate)
        end
        if existPT                                                          % power flux at top surface
            indT = Pz_top(:,i)>thresP;                                      % positive power flux at top
        end
    else                                                                    % use sign of vertical wavenumber
        if existkB                                                          % wavenumbers at bottom surface defined
            indB = all(squeeze(real(kz_bot(:,i,:)))<-thresK,2);             % negative wavenumbers at bottom (away from plate)
        end
        if existkT                                                          % wavenumbers at top surface
            indT = all(squeeze(real(kz_top(:,i,:)))>thresK,2);              % negative wavenumbers at top
        end
    end
    ind = indB & indT;
    indRemove = ~ind;
    
    % indRemove = indRemoveB | indRemoveT;                                    % combine both criteria

    k(indRemove,i) = nan;                                                   % remove unwanted solutions
    cp(indRemove,i) = nan;
    att(indRemove,i) = nan;
    if existkB
        kz_bot(indRemove,i,:) = nan;
    end
    if existkT
        kz_top(indRemove,i,:) = nan;
    end
    if existPB
        Pz_bot(indRemove,i) = nan;
    end
    if existPT
        Pz_top(indRemove,i) = nan;
    end
    if existCg
        cg(indRemove,i) = nan;
    end
    if existPhi
        phi(:,indRemove,i) = nan;
    end
    if existU
        u(:,indRemove,i) = nan;
    end
end

obj.kz_bottom = kz_bot;
obj.kz_top    = kz_top;
obj.Pz_bottom = Pz_bot;
obj.Pz_top    = Pz_top;
obj.k         = k;
obj.cp        = cp;
obj.att       = att;
obj.cg        = cg;
obj.phi       = phi;
obj.u         = u;

end
