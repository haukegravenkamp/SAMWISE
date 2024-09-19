function Delta = multipar_delta(A,indices)

%MULTIPAR_DELTA   Delta matrices for a multiparameter eigenvalue problem
%
% Delta = MULTIPAR_DELTA(A) returns set of k+1 Delta matrices, which are 
% operator determinants related to the multiparameter eigenvalue problem
%
% A{1,1} x1 = lambda(1) A{1,2} x1 +  ... + lambda(k) A{1,k+1} x1 
% A{2,1} x2 = lambda(1) A{2,2} x2 +  ... + lambda(k) A{2,k+1} x2 
% ...
% A{k,1} xk = lambda(1) A{k,2} xk +  ... + lambda(k) A{k,k+1} xk 
%
% See also: TWOPAR_DELTA, THREEPAR_DELTA
%
% If you want to compute just some of the matrices Delta, use 
% Delta = MULTIPAR_DELTA(A,indices), where indices is a vector with indices
% 0 to k of Delta matrices that you want to compute.

% MultiParEig toolbox
% B. Plestenjak and A. Muhic, University of Ljubljana
% FreeBSD License, see LICENSE.txt

% Last revision: 3.7.2022

    k = length(A) - 1; % number of parameters 
    if nargin<2
        indices = 0:k;
    end
    
    % computation of Delta matrices
    Delta = cell(1,k+1);
    for j = 1:length(indices)
        r = indices(j);
        Delta{r+1} = krondet(A, r);
    end

end

%------------------------------------------------------------------------

function krondelta =  krondet(A,index)
% recursive computation of Delta matrices - operator determinants

    if nargin == 2
        % we replace index-th column with the first (if index>0)
        n = length(A);
        if index==0
            ind = 2:n;
        else
            ind = [2:index 1 (index+2):n];
        end
        A = A(:, ind);
    end
    n = length(A);
    if n == 1
        krondelta = A{1,1}; 
        return;
    end
    deter = 0; 
    sgn = 1;
    for k = 1:n
        indj = [1:(k-1), (k+1):n];
        indi = 2:n;   
        deter = deter + sgn*kron( A{1, k}, krondet(A(indi, indj)) );
        sgn = -sgn;
    end
    
    krondelta = deter;
    
end
