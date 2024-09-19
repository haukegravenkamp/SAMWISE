
function shp = shapeFunctionsXnY(p,eta12,nModes,nDom,intTab,iT)

pDom = p(1:2);
pSid = p(3:end);

shp = zeros(3,nModes);                                                      % 2D shape functions

psiEta = @(eta) shpLagrange(2,eta,intTab);                                  % linear shape functions

eta1i = eta12(1);
eta2i = eta12(2);


indE = indicesCornersEdges(p,4);


%% COMPUTE SHAPE FUNCTIONS ALONG THE FOUR SIDES

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% collection of all shape functions along the sides
shpSides  = zeros(4,nModes);                                                % shape functions
shpSidesD = zeros(4,nModes);                                                % derivatives

etaSide = getEtaSide(eta1i,eta2i);                                          % relevant coordinates along the four edges
signDer = [1,1,-1,-1];                                                      % account for inverse direction of the sides 3 and 4

for iSide = 1:4                                                             % loop sides

    shpT = shpLagrange(pSid(iSide)+1,etaSide(iSide),intTab);
    shpSides(iSide,nonnan(indE(iSide,:)))  = shpT(1,:);
    shpSidesD(iSide,nonnan(indE(iSide,:))) = shpT(2,:)*signDer(iSide);
end


%% COMPUTE SHAPE FUNCTIONS AT THE FOUR CORNERS

% collection of all shape functions at the corners
shpCorners   = zeros(4,nModes);                                             % shape functions
shpCorners(:,1:4) = eye(4);

%% shape functions
psi1 = psiEta(eta1i);                                                       % projector in first direction
Peta1 = psi1(1,1)*shpSides(4,:) + psi1(1,2)*shpSides(2,:);

psi2 = psiEta(eta2i);                                                       % projector in second direction
Peta2 = psi2(1,1)*shpSides(1,:) + psi2(1,2)*shpSides(3,:);

psiProd = psi1(1,:)'*psi2(1,:);                                             % all products of blending functions
psiProd = psiProd([1 2 4 3]);

Peta12 = sum(psiProd.*shpCorners',2)';                                      % 2D interpolation of corner values
shp(1,:) = Peta1+Peta2-Peta12;                                              % shape functions


%% Derivatives, first direction

Peta1 = psi1(2,1)*shpSides(4,:) + psi1(2,2)*shpSides(2,:);                  % projector in first direction

Peta2 = psi2(1,1)*shpSidesD(1,:) + psi2(1,2)*shpSidesD(3,:);                % projector in second direction

psiProd = psi1(2,:)'*psi2(1,:);                                             % all products of blending functions
psiProd = psiProd([1 2 4 3]);

Peta12 = sum(psiProd.*shpCorners',2)';                                      % 2D interpolation of corner values

shp(2,:) = Peta1+Peta2-Peta12;                                              % shape functions

%% Derivatives, second direction

Peta1 = psi1(1,1)*shpSidesD(4,:) + psi1(1,2)*shpSidesD(2,:);                  % projector in first direction

Peta2 = psi2(2,1)*shpSides(1,:) + psi2(2,2)*shpSides(3,:);                    % projector in second direction

psiProd = psi1(1,:)'*psi2(2,:);                                               % all products of blending functions
psiProd = psiProd([1 2 4 3]);

Peta12 = sum(psiProd.*shpCorners',2)';                                        % 2D interpolation of corner values

shp(3,:) = Peta1+Peta2-Peta12;                                                % shape functions

%% bubble functions

shp1 = shpLagrange(pDom(1)+1,eta1i,intTab);
shp2 = shpLagrange(pDom(2)+1,eta2i,intTab);

% domain modes
ind = (size(shp,2)-nDom) + (1:nDom);
[N1,N2]= meshgrid(shp1(1,2:end-1),shp2(1,2:end-1));
shp(1,ind)=N1(:).*N2(:);
[N1,N2]= meshgrid(shp1(2,2:end-1),shp2(1,2:end-1));
shp(2,ind)=N1(:).*N2(:);
[N1,N2]= meshgrid(shp1(1,2:end-1),shp2(2,2:end-1));
shp(3,ind)=N1(:).*N2(:);

%% Apply blending to create nodal shape functions
if ~isempty(iT)
    shp(1,:) = shp(1,:)*iT;
    shp(2,:) = shp(2,:)*iT;
    shp(3,:) = shp(3,:)*iT;
end

    function eS = getEtaSide(e1,e2)
        eS = [e1,e2,-e1,-e2];                                % relevant coordinates along the four edges

    end



end
