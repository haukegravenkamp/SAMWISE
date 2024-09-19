function [k, evAll, X] = eig_Leaky_approximate(E0,E1,E2,M,R,typeCoupling,kappa,w,cAll)


kFree = sort(w./cAll);
kFree = [0;kFree;kFree(end)*1.1;kFree(end)*5];

k0 = (kFree(2:end)+kFree(1:end-1))/2;

nK = numel(k0);
kT{nK} = [];
kyT{nK} = [];
nEv = size(E1,1)*2;

nC = numel(kappa);

s = getSigns(nC);

for i = 1:nK
    kTest = [kFree(i),(kFree(i)+kFree(i+1))/2,kFree(i+1)];
    coeff = getInterpolation(kappa,kTest,typeCoupling);

    k = zeros(nEv*numel(s),1);
    ky = zeros(nEv*numel(s),nC);
    for iS = 1:size(s,1)
        [A0,A1,A2] = getRterms(R,coeff,s(iS,:),typeCoupling);
        [X1,k1]=polyeig(-E2+w^2*M + A0, 1i*E1 + A1, -E0 + A2);
        k((1:nEv) + nEv*(iS-1)) = k1;
        ky((1:nEv) + nEv*(iS-1),:) = s(iS,:).*sqrt(kappa.^2-k1.^2);
    end
    indR = (real(k)<0);
    indR = indR | (abs(real(k))) < kFree(i)*1  | (abs(real(k))) > kFree(i+1)*1.0;
    k(indR) = [];
    ky(indR,:) = [];
    kT{i} = k;
    kyT{i} = ky;

end
k = vertcat(kT{:});
evAll = [k,vertcat(kyT{:})];

    function coeff = getInterpolation(kappa,kTest,cType)
        switch cType
            case 'F'
                m = 1;
            case 'FF'
                m = [1;1];
            case 'S'
                m = repmat(kTest,2,1);
            case 'FS'
                m = [kTest;kTest;[1 1 1]];
            case 'SS'
                m = repmat(kTest,4,1);
        end

        coeff = zeros(numel(kappa),3);
        for iK = 1:numel(kappa)
            coeff(iK,:) = LagrangeQuadratic(kTest,m(iK,:).*sqrt(kappa(iK).^2-kTest.^2));
        end
    end

    function s = getSigns(nC)
        s = zeros(2^nC,nC);
        for j = 1:nC
            m = ones(2^(j-1),1);
            s(:,j) = repmat([m;-m],size(s,1)/numel(m)/2,1);
        end
    end
    function [A0,A1,A2] = getRterms(R,coeff,s,cType)
        switch cType
            case 'F'
                R = R{1};
                A0 = s*coeff(3)*1i*R;
                A1 = s*coeff(2)*1i*R;
                A2 = s*coeff(1)*1i*R;
            case 'S'
                A0 = s(1)*coeff(1,3)*R{1} + s(2)*coeff(2,3)*R{2};
                A1 = s(1)*coeff(1,2)*R{1} + s(2)*coeff(2,2)*R{2};
                A2 = s(1)*coeff(1,1)*R{1} + s(2)*coeff(2,1)*R{2};
            case 'FF'
                A0 = s(1)*coeff(1,3)*1i*R{1} + s(2)*coeff(2,3)*1i*R{2};
                A1 = s(1)*coeff(1,2)*1i*R{1} + s(2)*coeff(2,2)*1i*R{2};
                A2 = s(1)*coeff(1,1)*1i*R{1} + s(2)*coeff(2,1)*1i*R{2};
            case 'FS'
                A0 = s(1)*coeff(1,3)*R{1} + s(2)*coeff(2,3)*R{2} + s(3)*coeff(3,3)*1i*R{3};
                A1 = s(1)*coeff(1,2)*R{1} + s(2)*coeff(2,2)*R{2} + s(3)*coeff(3,2)*1i*R{3};
                A2 = s(1)*coeff(1,1)*R{1} + s(2)*coeff(2,1)*R{2} + s(3)*coeff(3,1)*1i*R{3};
            case 'SS'
                A0 = s(1)*coeff(1,3)*R{1} + s(2)*coeff(2,3)*R{2} + s(3)*coeff(3,3)*R{3} + s(4)*coeff(4,3)*R{4};
                A1 = s(1)*coeff(1,2)*R{1} + s(2)*coeff(2,2)*R{2} + s(3)*coeff(3,2)*R{3} + s(4)*coeff(4,2)*R{4};
                A2 = s(1)*coeff(1,1)*R{1} + s(2)*coeff(2,1)*R{2} + s(3)*coeff(3,1)*R{3} + s(4)*coeff(4,1)*R{4};

        end
    end

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
