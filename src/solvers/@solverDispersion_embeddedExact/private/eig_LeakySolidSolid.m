function [eigk,sols,X] = eig_LeakySolidSolid(L2,L1,L0,M,R,kL1,kS1,kL2,kS2,mu,opts)

% Returns eigenvalues and eigenvectors of the parameterized nonlinear eigenvalue problem
%
% [(ik)^2*L2 + ik*L1 + L0 + w^2 M + ik*ieta1*R11 + ik*ieta2*R21 + ik*ieta3*R12 + ik*ieta4*R22
%           + (ikL1)^2*R31 + (ikS1)^2*R41 + (ikL2)^2*R32 + (ikS2)*R42 ]*u = 0,
%
% where eta1 = sqrt(kL1^2-k^2), eta2 = sqrt(kS1^2-k^2), eta3 = sqrt(kS2^2-k^2), eta3 = sqrt(kL2^2-k^2),
% for a fixed mu, kL1, kS1, kL2 and kS2.
%
% We use nonlinearization and write the problem as a 6-parameter eigenvalue problem
%
% Ai1 x_i = (lambda*Ai2 + xi1*Ai3 + xi2*Ai4 + xi3*Ai5 + xi4*Ai6 + xi5*Ai7)x_i for i = 1,..,6
%
% where lambda = 1i*k, xi1 = 1i*k*1i*sqrt(kL1^2-k^2), xi2 = 1i*k*1i*sqrt(kS1^2-k^2), xi3 = (1i*k)^2,
% xi4 = 1i*k*1i*sqrt(kL2^2-k^2), xi5 = 1i*k*1i*sqrt(kS2^2-k^2),
%
% Output
%  - eigk: eigenvalues k
%  - sols: values [k eta1 eta2 eta3 eta4 k^2]

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% We use real computation for real matrices and parameters

% References
% 1) E. Ringh, E. Jarlebring: Nonlinearizing two-parameter eigenvalue
%    problems, SIAM J. Matrix Anal. Appl. 42 (2021) 775-799
% 2) B. Plestenjak: MultiParEig, MATLAB Central File Exchange.

% Bor Plestenjak 2023

if nargin<11
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
if ~isfield(opts,'refine'),    opts.refine = 0;  end
if ~isfield(opts,'solver'),    opts.solver = 'eig';  end
if ~isfield(opts,'twosideRQ'), opts.twosideRQ = 0;   end
if isfield(opts,'showgap'),    showgap   = opts.showgap;      else, showgap  = 0;  end
if isfield(opts,'maxlogres'),  maxlogres = opts.maxlogres;    else, maxlogres = 1e-3;  end

R11 = R{1};
R21 = R{2};
R12 = R{3};
R22 = R{4};

% A{1,1} = -(L0 + mu^2*M - kL1^2*R31 -kS1^2*R41  - kL2^2*R32 - kS2^2*R42);
A{1,1} = -(L0 + mu^2*M);
A{1,2} =  L1;
A{1,3} =  R11;
A{1,4} =  R21;
A{1,5} =  L2;
A{1,6} =  R12;
A{1,7} =  R22;

% 2 x 2 equation introduces xi1^2 = -xi3*(kL1^2+xi3), where xi3 = lambda^2
A{2,1} = -[0 -kL1^2;0 0];
A{2,2} = zeros(2);
A{2,3} = eye(2);
A{2,4} = zeros(2);
A{2,5} = [0 -1;1 0];
A{2,6} = zeros(2);
A{2,7} = zeros(2);

% 2 x 2 equation introduces xi2^2 = -xi3*(kS1^2+xi3), where xi3 = lambda^2
A{3,1} = -[0 -kS1^2;0 0];
A{3,2} = zeros(2);
A{3,3} = zeros(2);
A{3,4} = eye(2);
A{3,5} = [0 -1;1 0];
A{3,6} = zeros(2);
A{3,7} = zeros(2);

% 2 x 2 equation introduces xi3^2 = lambda^2
A{4,1} = -[0 0;0 1];
A{4,2} = [0 1; 1 0];
A{4,3} = zeros(2);
A{4,4} = zeros(2);
A{4,5} = [1 0; 0 0];
A{4,6} = zeros(2);
A{4,7} = zeros(2);

% 2 x 2 equation introduces xi4^2 = -xi3*(kL2^2+xi3), where xi3 = k^2
A{5,1} = -[0 -kL2^2;0 0];
A{5,2} = zeros(2);
A{5,3} = zeros(2);
A{5,4} = zeros(2);
A{5,5} = [0 -1;1 0];
A{5,6} = eye(2);
A{5,7} = zeros(2);

% 2 x 2 equation introduces xi5^2 = -xi3*(kS2^2+xi3), where xi3 = k^2
A{6,1} = -[0 -kS2^2;0 0];
A{6,2} = zeros(2);
A{6,3} = zeros(2);
A{6,4} = zeros(2);
A{6,5} = [0 -1;1 0];
A{6,6} = zeros(2);
A{6,7} = eye(2);

if use_shift
    % we shift matrices A{j,5} to make Delta_0 nonsingular
    for j =1:6
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
        % tmp_lambda = tmp_lambda(ord(infeigs+1:end),:);
        ind = ord(infeigs+1:end);
        tmp_lambda = tmp_lambda(ind,:);
        if computeVectors
            X = X(ind,:);
        end
        if showgap
            fprintf('Keep %d solutions, dif: %5.5e %5.5e \n',length(tmp_lambda),ll(infeigs),ll(infeigs+1))
        end
    end

    % correct eigenvalues due to the shift
    full_lambda = tmp_lambda;
    for j = 1:6
        full_lambda(:,j) = tmp_lambda(:,j)./(1-shift*tmp_lambda(:,4));
    end
else
    full_lambda = tmp_lambda;
end

eigk = full_lambda(:,1)/1i;
sols = [eigk full_lambda(:,2)./eigk full_lambda(:,3)./eigk full_lambda(:,5)./eigk full_lambda(:,6)./eigk full_lambda(:,4)];
if computeVectors
    X = X(:,[1 2 3 5 6 4]);
end

% We take only solutions where computed eta1, eta2, eta3, eta4 give small residual
if maxlogres>0
    eta1_er =  abs(log(abs(sqrt(kL1^2-eigk.^2)./sols(:,2))));
    eta2_er =  abs(log(abs(sqrt(kS1^2-eigk.^2)./sols(:,3))));
    eta3_er =  abs(log(abs(sqrt(kL2^2-eigk.^2)./sols(:,4))));
    eta4_er =  abs(log(abs(sqrt(kS2^2-eigk.^2)./sols(:,5))));
    test_sr = eta1_er + eta2_er + eta3_er + eta4_er;
    tmpres = sort(test_sr);
    ind = find(test_sr<maxlogres);
    select = length(ind);
    if (select < length(eigk)) && showgap  && select >0
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
