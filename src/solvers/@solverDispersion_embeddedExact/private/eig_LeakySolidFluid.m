function [eigk,sols,X] = eig_LeakySolidFluid(L2,L1,L0,M,Rp,Rs,Rf,kL1,kS1,kF,mu,opts)

% Returns eigenvalues and eigenvectors of the parameterized nonlinear eigenvalue problem
%
% [(ik)^2*L2 + ik*L1 + L0 + w^2 M + ik*ibeta*R11 + ik*ieta*R21 + igamma*RF + (ikL1)^2*R31 + (ikS1)^2*R41 ]*u = 0,
%
% where beta = sqrt(kL1^2-k^2), eta = sqrt(kS1^2-k^2), for a fixed mu, kL1, kS1 and kF
%
% We use nonlinearization and write the problem as a 5-parameter eigenvalue problem
%
% Ai1*x_i = (lambda*Ai2 + xi1*Ai3 + xi2*Ai4 + xi3*Ai5 + xi4*Ai6)*x_i for i = 1,..,5
%
% where lambda = ik, xi1 = k*sqrt(kL1^2-k^2), xi2 = k*sqrt(kS1^2-k^2), xi3 = k^2, xi4 = sqrt(kF^2-k^2)
%
% Output
%  - eigk: eigenvalues k
%  - sols: values [k beta eta gamma]

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
% References
% 1) E. Ringh, E. Jarlebring: Nonlinearizing two-parameter eigenvalue
%    problems, SIAM J. Matrix Anal. Appl. 42 (2021) 775-799
% 2) B. Plestenjak: MultiParEig, MATLAB Central File Exchange.

% Bor Plestenjak 2023

if nargin<12
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

% A{4,1} = -(L0 + mu^2*M - k1^2*R3 -k2^2*R4);
A{5,1} = -(L0 + mu^2*M);
A{5,2} =  L1;
A{5,3} =  Rp;
A{5,4} =  Rs;
A{5,5} =  Rf;
A{5,6} =  L2;

% 2 x 2 equation introduces xi1^2 = -xi3*(k1^2+xi3), where xi3 = lambda^2
A{2,1} = -[0 -kL1^2 ;0 0];
A{2,2} = zeros(2);
A{2,3} = eye(2);
A{2,4} = zeros(2);
A{2,5} = zeros(2);
A{2,6} = [0 -1;1 0];

% 2 x 2 equation introduces xi2^2 = -xi3*(k2^2+xi3), where xi3 = lambda^2
A{3,1} = -[0 -kS1^2;0 0];
A{3,2} = zeros(2);
A{3,3} = zeros(2);
A{3,4} = eye(2);
A{3,5} = zeros(2);
A{3,6} = [0 -1;1 0];

A{4,1} = -[0 -kF^2; 1 0];
A{4,2} = zeros(2);
A{4,3} = zeros(2);
A{4,4} = zeros(2);
A{4,5} = [1 0; 0 1];
A{4,6} = [0 -1; 0 0];

% 2 x 2 equation introduces xi3^2 = lambda^2
A{1,1} = -[0 0;0 1];
A{1,2} = [0 1; 1 0];
A{1,3} = zeros(2);
A{1,4} = zeros(2);
A{1,5} = zeros(2);
A{1,6} = [1 0; 0 0];





if use_shift
    % we shift matrices A{j,6} to make Delta_0 nonsingular
    for j =1:5
        A{j,6} = A{j,6} + shift*A{j,1};
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
        [ll,ord] = sort(abs(1-shift*tmp_lambda(:,5)));
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
    for j = 1:5
        full_lambda(:,j) = tmp_lambda(:,j)./(1-shift*tmp_lambda(:,5));
    end
else
    full_lambda = tmp_lambda;
end

% Eigenvalues in full_lambda are of the form [ik, -k*beta, -k*eta, i*gamma, (ik)^2+kF^2]

% The output is [k, beta, eta, gamma]
eigk = full_lambda(:,1)/1i;
sols = [eigk full_lambda(:,2)./eigk full_lambda(:,3)./eigk full_lambda(:,4)/1i];
if computeVectors
    X = X(:,[5 2 3 4 1]);
end

% We take only solutions where computed beta, eta, gamma give small residual
if maxlogres>0
    beta_er =  abs(log(abs(sqrt(kL1^2-eigk.^2)./sols(:,2))));
    eta_er =   abs(log(abs(sqrt(kS1^2-eigk.^2)./sols(:,3))));
    gamma_er = abs(log(abs(sqrt(kF^2-eigk.^2)./sols(:,4))));

    test_sr = beta_er + eta_er + gamma_er;
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
