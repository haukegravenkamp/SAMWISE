function [proriolOrders,iV] = prepareInnerModesTri(p,intTab)

%% Prepare computation of inner modes

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% NODAL POSITIONS ACCORDING TO POZRIKIDIS TEXTBOOK, SECTION 5.8
etaGLL = intTab(p+1).GLL.xi;                                                     % nodes at GLL points
v = (etaGLL+1)/2;                                                                % scale to 0->1

% LOCAL NODE COORDINATES ON TRIANGLE (BETWEEN 0 AND 1!)
xiCor  = [0; 1; 0];                                                         % corner nodes
etaCor = [0; 0; 1];

xiEdge  = [v(2:end-1),v(end-1:-1:2),zeros(1,numel(v)-2)]';                  % edge nodes
etaEdge = [zeros(1,numel(v)-2),v(2:end-1),v(end-1:-1:2)]';

[xiDom,etaDom] = triIntCoordsLocal(p,etaGLL);

xi  = [xiCor;  xiEdge;  xiDom ];                                            % collect all modes in correct order
eta = [etaCor; etaEdge; etaDom];

xid  = 2*xi./(1-eta)-1;                                                     % map on standard square ("xi dash", "eta dash")
etad = 2*eta-1;

xid(isnan(xid)) = 0;                                                        % For eta=1 the transformation to xi is undefined. However, the proriol

nnodes = sum(1:(p+1));                                                      % number of nodes
V = zeros(nnodes);                                                          % allocate Vandermonde matrix

proriolOrders = zeros(nnodes,2);                                            % orders of the required proriol polynomials
c = 1;                                                                      % counter
for k = 0:p
    for l = 0:p-k
        proriolOrders(c,:) = [k,l];                                         % orders of current proriol polynomial
        V(c,:) = proriol(k,l,xid,etad)';                                    % compute one row of Vandermonde matrix
        c = c+1;                                                            % update counter
    end
end
iV = V\eye(size(V));                                                        % inverse of Vandermonde matrix

end


