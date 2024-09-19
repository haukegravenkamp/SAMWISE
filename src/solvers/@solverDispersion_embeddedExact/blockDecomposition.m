function [coeffM,R,dofD] = blockDecomposition(obj,coeffM,R,typeC,opt)

E0  = coeffM{1};
E11 = coeffM{2};
E12 = coeffM{3};
E2  = coeffM{4};
C0  = coeffM{5};
M0  = coeffM{6};

decomposeEVP = opt.numerics.decomposeEVP;
if ~(strcmp(typeC,'F') || strcmp(typeC,'FF'))
    decomposeEVP = false;
end

if decomposeEVP

    if strcmp(typeC,'F') && ~obj.sameUnbMedia
        indPl = 1:size(E2,1)-1;
        iU    = size(E2,1);
    elseif strcmp(typeC,'FF') || strcmp(typeC,'F')
        indPl = 1:size(E2,1)-2;
        iU    = size(E2,1)-2 + (1:2);
    end
    [ind,nBl] = solverDispersion.blockDecomposition(E2(indPl,indPl)+M0(indPl,indPl),1e-6);
    E1 = E12-E11;
    AC = E0 + E2 + M0 + C0;
end
dofD = false;
if decomposeEVP && (nBl == 2) ...                                           % two blocks (plane strain or plane stress)
        && ~any(E1([ind{1}],[ind{1}]),'all') ...
        && ~any(E1([ind{2}],[ind{2}]),'all') ...
        && ~any(AC([ind{1}],[ind{2}]),'all')...
        && ~any(AC([ind{2}],[ind{1}]),'all')

    i1 = ind{1};
    i2 = ind{2};
    iA = [i1,i2,iU];

    M0 = M0(iA,iA);
    for iR = 1:numel(R)
        R{iR} = R{iR}(iA,iA);
    end

    E0 = [E0(i1,iA);  E1(i2,i1), E0(i2,[i2,iU]);E0(iU,iA);];
    E2 = [E2(i1,i1), -E1(i1,[i2,iU]); E2([i2,iU],iA)];
    E11 = zeros(size(E0));
    E12 = zeros(size(E0));

    [~,dofD]=sort(iA);                                                      % indices to revert decomposition
    dofD(2,end) = 0;                                                        % store indices of dofs multiplied by ik in the same array
    dofD(2,i2) = 1;
    dofD(2,iU) = 1;

elseif decomposeEVP && (nBl == 3) ...                                       % three blocks (including out-of-plane modes)
        && ~any(E1([ind{1}],[ind{1}]),'all') ...
        && ~any(E1([ind{2}],[ind{2}]),'all') ...
        && ~any(E1([ind{3}],[ind{3}]),'all') ...
        && ~any(AC([ind{1}],[ind{2}]),'all')...
        && ~any(AC([ind{1}],[ind{3}]),'all')...
        && ~any(AC([ind{2}],[ind{1}]),'all')...
        && ~any(AC([ind{2}],[ind{3}]),'all')...
        && ~any(AC([ind{3}],[ind{1}]),'all')...
        && ~any(AC([ind{3}],[ind{2}]),'all')

    if ~any(E1([ind{1}],[ind{:}]),'all')                                    % find the isolated block of out-of-plane modes and move it to the end for consistency
        i1 = ind{2};
        i2 = ind{3};
        i3 = ind{1};
    elseif ~any(E1([ind{2}],[ind{:}]),'all')
        i1 = ind{1};
        i2 = ind{3};
        i3 = ind{2};
    elseif ~any(E1([ind{3}],[ind{:}]),'all')
        i1 = ind{1};
        i2 = ind{2};
        i3 = ind{3};
    end

    iA = [i1,i2,i3,iU];


    M0 = M0(iA,iA);
    for iR = 1:numel(R)
        R{iR} = R{iR}(iA,iA);
    end

    E0 = [E0(i1,iA);  E1(i2,i1), E0(i2,[i2,i3,iU]);E0([i3,iU],iA);];
    E2 = [E2(i1,i1), -E1(i1,[i2,i3,iU]); E2([i2,i3,iU],iA)];
    E11 = zeros(size(E0));
    E12 = zeros(size(E0));

    [~,dofD]=sort(iA);                                                      % indices to revert decomposition
    dofD(2,end) = 0;                                                        % store indices of dofs multiplied by ik in the same array
    dofD(2,i2) = 1;
    dofD(2,iU) = 1;

end

coeffM = {E0,E11,E12,E2,C0,M0};

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
