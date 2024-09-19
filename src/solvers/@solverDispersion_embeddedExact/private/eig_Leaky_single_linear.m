function [eigk,full_eig,X] = eig_Leaky_single_linear(L2,L0,M,R,k1,mu,dofD,opts)
% Similar to eig_Leaky_single but without linear term in k -> no L1 matrix
% 
% Returns eigenvalues and eigenvectors of the parameterized nonlinear eigenvalue problem
% 
% [(1i*k)^2*L2 + L0 + w^2 M + 1i*beta*R]*u = 0,
%
% where beta = sqrt(k1^2-k^2), for a fixed mu and k1.
%
% We use nonlinearization and write the problem as a 2-parameter eigenvalue problem
%
% Ai1 x_i = (lambda*Ai2 + xi1*Ai3 + xi2*Ai4)x_i for i = 1,2
% 
% where lambda = 1i*k, xi1 = 1i*sqrt(k1^2-k^2), xi2 = lambda^2 + k1^2
%
% Output
%  - eigk: eigenvalues k
%  - full_eig: values [k beta] 

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% References
% 1) E. Ringh, E. Jarlebring: Nonlinearizing two-parameter eigenvalue
%    problems, SIAM J. Matrix Anal. Appl. 42 (2021) 775-799
% 2) B. Plestenjak: MultiParEig, MATLAB Central File Exchange. 

% Bor Plestenjak 2023

if nargin<8
    opts = [];
end
if nargout>2
    computeVectors = true;
else
    computeVectors = false;
end

if isfield(opts,'use_shift'),  use_shift = opts.use_shift;    else, use_shift = 1;  end
if isfield(opts,'shift'),      shift     = opts.shift;        else, shift     = 1;  end
if isfield(opts,'infeigs'),    infeigs   = opts.infeigs;      else, infeigs   = 4;  end
if ~isfield(opts,'rand_orth'), opts.rand_orth = 1; end
if ~isfield(opts,'refine'),    opts.refine = 0;  end
if ~isfield(opts,'solver'),    opts.solver = 'eig';  end
if ~isfield(opts,'twosideRQ'), opts.twosideRQ = 0;  end
if isfield(opts,'showgap'),    showgap   = opts.showgap;      else, showgap  = 0;  end

% If all matrices and parameters are real, we get a real MEP that is solved faster

% this 2 x 2 equation gives xi1^2 = -(lambda^2 + k1^2) = - xi3 - k1^2
A{1,1} = -[0 -k1^2; 1 0];
A{1,2} = [1 0; 0 1];
A{1,3} = [0 -1; 0 0];

% original equation
A{2,1} = -(L0 + mu^2*M);
A{2,2} = R;
A{2,3} = L2;

% shift is not so accurate as singular, but it is faster
if use_shift 
    % we shift matrices A{j,4} to make Delta_0 nonsingular
    for j = 1:2
        A{j,3} = A{j,3} + shift*A{j,1};
    end
else
    % solve a singular MEP
    opts.singular = 1;
end

if ~computeVectors
   tmp_lambda = multipareig(A,opts);
else
   [tmp_lambda,X] = multipareig(A,opts);
end

if use_shift
    if infeigs>0
        % we discard infeigs eigenvalues with the smallest value of 
        % 1-shift*xi3, which should be zero or close to zero for infinite eigenvalues
        [ll,ord] = sort(abs(1-shift*tmp_lambda(:,2)));
        ind = ord(infeigs+1:end);
        tmp_lambda = tmp_lambda(ind,:);
        if computeVectors
        X = X(ind,:);
        end
        if showgap
            fprintf('Keep %d solutions, dif: %5.5e %5.5e \n',length(tmp_lambda),ll(infeigs),ll(infeigs+1))
        end
    end

    % transform the eigenvalues due to the shift
    full_eig = tmp_lambda;
    for j = 1:2
        full_eig(:,j) = tmp_lambda(:,j)./(1-shift*tmp_lambda(:,2));
    end
else
    full_eig = tmp_lambda;
end
                                       
kU = -1i*full_eig(:,1);                                                     % unbounded wavenumber
kU = repmat(kU,2,1);
eigk = sqrt(-full_eig(:,2));
eigk = [eigk;-eigk];
full_eig = [eigk,kU];

if computeVectors
    X = X(:,[2 1]);
    X = [X;X];
end

% correct eigenvectors in case of decomposition
i2 = logical(dofD(2,:));
iS = dofD(1,:);
for i = 1:numel(eigk)
    if computeVectors
        X{i,1} = X{i,1}(iS);
        X{i,1}(i2) =  -1i*X{i,1}(i2)/eigk(i);
    end
end

end                                            

