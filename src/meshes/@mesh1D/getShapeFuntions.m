function [N,Nd1,Nd2,Na,Nad1,Nad2] = getShapeFuntions(p,eta,dof,intTab,~,~,~)

nEta = numel(eta);

N = zeros(nEta,(p+1));                                                     % allocate scalar shape functions
Nd1 = N;
Nd2 = N;

Na = zeros(dof,(p+1)*dof,nEta);                                             % allocate vectorial shape functions for all dofs
Nad1 = Na;
Nad2 = Na;

for i = 1:nEta

    shp = shpLagrange(p+1,eta(i),intTab);                                   % get shape functions and derivatives

    N(i,:)   = shp(1,:);                                                    % shape function matrix
    Nd1(i,:) = shp(2,:);                                                    % shape function derivatives

    Na(:,:,i)   = kron(shp(1,:),eye(dof));                                  % shape function matrix for all dofs
    Nad1(:,:,i) = kron(shp(2,:),eye(dof));                                  % shape function derivatives

end
dofSort = [1:dof,p*dof+(1:dof),(dof+1):p*dof];
Na = Na(:,dofSort,:);
Nad1 = Nad1(:,dofSort,:);

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
