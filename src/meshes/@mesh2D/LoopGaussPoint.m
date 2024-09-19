function eleMat = LoopGaussPoint(N,Ndeta1,Ndeta2,bxL,byL,bzL,~,bxR,byR,bzR,~,coeff,eta,~,~,J,Jdet,wgts,nDofs,nM,matCompute)

eleMat = zeros(nDofs,nDofs,nM);                                             % allocate fe-matrices for current term (stacked in case of several matrices)

ydeta1 = J(2,2,:);
zdeta1 = J(2,3,:);
ydeta2 = J(3,2,:);
zdeta2 = J(3,3,:);

for iGP = 1:size(eta,1)                                                      % Gauss point loop



    wi = wgts(iGP);                                                         % current integration weight
    Jeti = Jdet(iGP);                                                       % Jacobian determinant at current Gauss point

    B11 = bxL*N(:,:,iGP);                                                   % B-matrices
    B12 = bxR*N(:,:,iGP);
    B21 = 1/Jeti*( zdeta2(iGP)*byL - ydeta2(iGP)*bzL)*Ndeta1(:,:,iGP) + ...
          1/Jeti*(-zdeta1(iGP)*byL + ydeta1(iGP)*bzL)*Ndeta2(:,:,iGP);
    B22 = 1/Jeti*( zdeta2(iGP)*byR - ydeta2(iGP)*bzR)*Ndeta1(:,:,iGP) + ...
          1/Jeti*(-zdeta1(iGP)*byR + ydeta1(iGP)*bzR)*Ndeta2(:,:,iGP);

    c = coeff(:,:,iGP);                                                     % coefficient (material parameter) at Gauss point

    if matCompute(1)
        eleMat(:,:,1) = eleMat(:,:,1) + wi*B11'*c*B12*Jeti;                 % first term (no spatial derivatives)
    end
    if matCompute(2)                                                        % first derivatives in trial functions
        eleMat(:,:,2) = eleMat(:,:,2) + wi*B21'*c*B12*Jeti;
    end
    if matCompute(3)                                                        % first derivatives in test functions
        eleMat(:,:,3) = eleMat(:,:,3) + wi*B11'*c*B22*Jeti;
        % requires modification should cases with derivatives in only test
        % or trial functions become relevant
    end
    if matCompute(4)                                                        % first derivatives in both
        eleMat(:,:,4) = eleMat(:,:,4) + wi*B21'*c*B22*Jeti;
    end
end






end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
