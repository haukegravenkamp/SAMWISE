function [eigk,sols,X] = eig_LeakySolid(L2,L1,L0,M,R1,R2,k1,k2,mu,opts)

% Returns eigenvalues and eigenvectors of the parameterized nonlinear eigenvalue problem
%
% [(1i*k)^2*L2 + 1i*k*L1 + L0 + w^2 M + 1i*k*1i*beta*R1 + 1i*k*1i*eta*R2 + (1i*k_1)^2*R3 + (1i*k_2)^2*R4 ]*u = 0,
%
% where beta = sqrt(k1^2-k^2), eta = sqrt(k2^2-k^2), for a fixed mu, k_1 and k_2.
%
% We use nonlinearization and write the problem as a 4-parameter eigenvalue problem
%
% Ai1 x_i = (lambda*Ai2 + xi1*Ai3 + xi2*Ai4 + xi3*Ai5)x_i for i = 1,..,4
%
% where lambda = 1i*k, xi1 = lambda*1i*sqrt(k1^2+lambda^2), xi2 = lambda*1i*sqrt(k2^2+lambda^2), xi3 = lambda^2
%
% Output
%  - eigk: eigenvalues k
%  - sols: values [k beta eta]

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% References
% 1) E. Ringh, E. Jarlebring: Nonlinearizing two-parameter eigenvalue
%    problems, SIAM J. Matrix Anal. Appl. 42 (2021) 775-799
% 2) B. Plestenjak: MultiParEig, MATLAB Central File Exchange.

% Bor Plestenjak 2023

if nargin<10
    opts = [];
end
if nargout>2
    computeVectors = true;
else
    computeVectors = false;
end

if isfield(opts,'use_shift'),  use_shift = opts.use_shift;    else, use_shift = 1;  end
if isfield(opts,'shift'),      shift     = opts.shift;        else, shift     = 1;  end
if isfield(opts,'infeigs'),    infeigs   = opts.infeigs;      else, infeigs   = 0;  end
if ~isfield(opts,'rand_orth'), opts.rand_orth = 1; end
if ~isfield(opts,'solver'),    opts.solver = 'eig';  end
if ~isfield(opts,'refine'),    opts.refine = 0;  end
if ~isfield(opts,'twosideRQ'), opts.twosideRQ = 1;  end
if isfield(opts,'showgap'),    showgap   = opts.showgap;      else, showgap  = 0;  end
if isfield(opts,'maxlogres'),  maxlogres = opts.maxlogres;    else, maxlogres = 1e-2;  end

% A{4,1} = -(L0 + mu^2*M - k1^2*R3 -k2^2*R4);
A{4,1} = -(L0 + mu^2*M);
A{4,2} =  L1;
A{4,3} =  R1;
A{4,4} =  R2;
A{4,5} =  L2;

% 2 x 2 equation introduces xi1^2 = -xi3*(k1^2+xi3), where xi3 = lambda^2
A{2,1} = -[0 -k1^2 ;0 0];
A{2,2} = zeros(2);
A{2,3} = eye(2);
A{2,4} = zeros(2);
A{2,5} = [0 -1;1 0];

% 2 x 2 equation introduces xi2^2 = -xi3*(k2^2+xi3), where xi3 = lambda^2
A{3,1} = -[0 -k2^2;0 0];
A{3,2} = zeros(2);
A{3,3} = zeros(2);
A{3,4} = eye(2);
A{3,5} = [0 -1;1 0];

% 2 x 2 equation introduces xi3^2 = lambda^2
A{1,1} = -[0 0;0 1];
A{1,2} = [0 1; 1 0];
A{1,3} = zeros(2);
A{1,4} = zeros(2);
A{1,5} = [1 0; 0 0];

if use_shift
    % we shift matrices A{j,5} to make Delta_0 nonsingular
    for j =1:4
        A{j,5} = A{j,5} + shift*A{j,1};
    end
else
    % problem is singular
    opts.singular = 1;
end

if ~computeVectors
    tmp_lambda = multipareig(A,opts);
else
    [tmp_lambda,X] = multipareig(A,opts);
end

if use_shift
    if infeigs>0
        % in this case we discard infeigs infinite eigenvalues, we take
        % infeigs with the smallest values of |1-shift*tmp_lambda(:,4)|,
        % because this should be 0 for infinite eigenvalues
        [ll,ord] = sort(abs(1-shift*tmp_lambda(:,4)));
        ind = ord(infeigs+1:end);
        tmp_lambda = tmp_lambda(ind,:);
        X = X(ind,:);
        if showgap
            fprintf('Keep %d solutions, dif: %5.5e %5.5e \n',length(tmp_lambda),ll(infeigs),ll(infeigs+1))
        end
    end

    % correct eigenvalues due to the shift
    full_lambda = tmp_lambda;
    for j = 1:4
        full_lambda(:,j) = tmp_lambda(:,j)./(1-shift*tmp_lambda(:,4));
    end
else
    full_lambda = tmp_lambda;
end

eigk = full_lambda(:,1)/1i;
sols = [eigk full_lambda(:,2)./eigk full_lambda(:,3)./eigk];

if computeVectors
    X = X(:,[4 2 3 1]);
end

% We take only solutions where beta and eta are close to sqrt(k1^2-k^2) and sqrt(k1^2-k^2)
if maxlogres>0
    test_sr = abs(log(abs(sqrt(k1^2-sols(:,1).^2)./sols(:,2)))) + abs(log(abs(sqrt(k2^2-sols(:,1).^2)./sols(:,3))));
    tmpres = sort(test_sr);
    ind = find(test_sr<maxlogres);
    select = length(ind);
    if (select < length(eigk)) && showgap && select >0
        fprintf('Keep %d out of %d solutions, gap: %5.5e %5.5e \n',select,length(tmp_lambda),tmpres(select),tmpres(select+1))
    end

    if ~isempty(ind)
        eigk = eigk(ind);
        sols = sols(ind,:);
        if computeVectors
            X = X(ind,:);
        end
    else
        eigk = nan;
        sols = sols(1,:)*nan;
        if computeVectors
            for i = 1:size(X,2)
                X{1,i}=X{1,i}*nan;
            end
            X=X(1,:);
        end
    end
end
