function obj = computePowerFlux(obj,bcd,mat)
PzB = zeros(size(obj.phi,2),size(obj.phi,3));                               % allocate vertical power flux
PzT = PzB;

existBT = [false,false];

for ib = 1:numel(bcd)                                                       % loop boundary conditions

    if isa(bcd(ib),'bcUnbounded')                                           % is a b.c. for unbounded domain
        matNoUnb = bcd(ib).material;                                        % material number of unbounded medium
        dofsPla = bcd(ib).dofs;                                             % dofs of plate
        dofsUnb = bcd(ib).dofsUnbounded;                                    % dofs of unbounded domain
        if isscalar(dofsPla)                                                % only 1 dof used of plate -> no vertical component
            continue
        end
        if strcmp(bcd(ib).location,'bottom')
            top = false;
            existBT(1) = true;
        elseif strcmp(bcd(ib).location,'top')
            top = true;
            existBT(2) = true;
        end

        if strcmp(mat(matNoUnb).behavior,'acoustic')                        % check whether coupling is to acoustic or elastic medium

            for iW = 1:size(obj.phi,3)
                w = obj.omega(iW);
                for iK = 1:size(obj.phi,2)

                    p = obj.phi(dofsUnb,iK,iW);
                    u = obj.phi(dofsPla,iK,iW);
                    v = -1i*w*u;
                    Poy = real(v'*p)/2;
                    if top
                        PzT(iK,iW) = Poy(2);
                    else
                        PzB(iK,iW) = Poy(2);
                    end
                end
            end

        elseif strcmp(mat(matNoUnb).behavior,'elastic')


            cl = mat(matNoUnb).elastic.cl;
            cs = mat(matNoUnb).elastic.cs;
            rMatrices = bcd(ib).rMatrices;
            r0 = rMatrices{1};
            r1 = rMatrices{2};
            r2 = rMatrices{3};
            r3 = rMatrices{4};
            r4 = rMatrices{5};

            if top
                kzAll = obj.kz_top;
            else
                kzAll = obj.kz_bottom;
            end

            for iW = 1:size(obj.phi,3)
                w = obj.omega(iW);
                kappaL = w/cl;
                kappaS = w/cs;
                for iK = 1:size(obj.phi,2)
                    c = obj.phi(dofsUnb,iK,iW);
                    u = obj.phi(dofsPla,iK,iW);
                    v = -1i*w*u;
                    k = obj.k(iK,iW);
                    ky = squeeze(kzAll(iK,iW,:));

                    tau = (k^2*r0 + k*ky(1)*r1 + k*ky(2)*r2 +...
                        kappaL^2*r3 + kappaS^2*r4)*c;

                     Poy = -real(v'*tau)/2;
                    
                     if top
                         PzT(iK,iW) = Poy;
                     else
                         PzB(iK,iW) = Poy;
                     end


                end
            end

        end

    end

end
if existBT(1)
    obj.Pz_bottom = PzB;
end
if existBT(2)
    obj.Pz_top = PzT;
end


end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
