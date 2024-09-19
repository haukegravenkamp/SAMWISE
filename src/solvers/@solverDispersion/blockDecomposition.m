%% decomposition of matrix flow E using eigenvectors Phi with accuracy thB
function [ind,nBl]=blockDecomposition(B,thB)
% apply transformation
normB = normest(B,1e-2);
B(abs(B)/normB<thB)=0;
% neglect small values
[p,~,r,~,~,~] = dmperm(B);
% permutation
nBl = numel(r)-1;
% number of blocks
ind=cellfun(@(i)p(r(i):r(i+1)-1),num2cell(1:nBl), 'UniformOutput' ,false);
% store block indices
end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
