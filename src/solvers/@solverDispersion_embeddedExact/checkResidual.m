function obj = checkResidual(obj,E0,E1,E2,M,typeC,RA,aDofs)


wA = obj.omega;
kA = obj.k;
phiA = obj.phi;
ky_top = obj.ky_top;
ky_bottom = obj.ky_bottom;

nF = numel(wA);
nM = size(kA,1);

res = nan(nM,nF);

for iF = 1:nF
    w = wA(iF);
    for iM = 1:nM
        k = kA(iM,iF);
        phi = phiA(aDofs,iM,iF);
        % kappa = sqrt(w^2./c.^2-k^2)
        if nF > 1
            kappaT = squeeze(ky_top(iM,iF,:));
        else
            kappaT = ky_top(iM,:);
        end
        if nF > 1
             kappaB = squeeze(ky_bottom(iM,iF,:));
        else
             kappaB = ky_bottom(iM,:);
        end
        if isempty(kappaB)
            kappa = kappaT;
        else
            kappa = kappaB;
        end
        switch typeC
            case 'S'                                                        % two solid unbounded domains
                R = 0;
                for i = 1:2
                    R = R + kappa(i)*k*RA{i};
                end

            case 'SS'                                                       % two solid unbounded domains
                R = 0;
                for i = 1:2                                                 % bottom
                    R = R + kappaB(i)*k*RA{i};
                end
                for i = 1:2                                                 % top
                    R = R + kappaT(i)*k*RA{i+2};
                end

            case 'F'                                                        % one fluid unbounded domain
                R =  kappa*1i*RA{1};

            case 'FF'                                                       % two fluid unbounded domains
                R =      kappaB*1i*RA{1};
                R =  R + kappaT*1i*RA{2};

            case 'FS'
                if isscalar(kappaB)                                         % bottom fluid
                    kappa = [kappaT;kappaB];                                % always solid first
                else
                    kappa = [kappaB;kappaT];
                end
                R = 0;
                for i = 1:2                                                 % bottom
                    R = R + kappa(i)*k*RA{i};
                end
                R =  R + kappa(3)*1i*RA{3};                                 % CHECK THIS SIGN


        end

        res(iM,iF) = norm((-k^2*E0 + 1i*k*E1 - E2 + w^2*M + R)*phi);

    end

end

textStr = 'min/max/mean residual of all solutions: ';
nstr1 = num2str(round(min(nonnan(res(:))),2,'significant'),'%.1d');
nstr2 = num2str(round(max(nonnan(res(:))),2,'significant'),'%.1d');
nstr3 = num2str(round(mean(nonnan(res(:))),2,'significant'),'%.1d');

disp([textStr, nstr1,', ', nstr2,', ', nstr3])

obj.res = res;


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
