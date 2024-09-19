
function plotShapefunctions(obj,pAll,nV,nDom)

if nargin < 2
    pAll = [1 1 1 1 1 1];
end
intTab = inttable(max(pAll));

N = 50;

eta1 = linspace(-1,1,N);
eta2 = linspace(-1,1,N);

[eta1, eta2] = meshgrid(eta1,eta2);
eta = [eta1(:),eta2(:)];

dof = 1;

% nDom = (pAll(1)-1)*(pAll(2)-1);
nEleNodes =  sum(pAll(3:end))+nDom;

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
[N,Ndeta1,Ndeta2] = obj.getShapeFuntions(pAll,eta,dof,intTab,nEleNodes,nDom,nV);

if nV == 3
    u=zeros(size(eta1));
    v=zeros(size(eta1));
    w=zeros(size(eta1));
    y=zeros(size(eta1));
    z=zeros(size(eta1));
    vert = [0 0; 1 0; 0 1];
    for i = 1:numel(eta1)
    [u(i),v(i),w(i)]=local2Barycentric(eta1(i),eta2(i));
    [y(i),z(i)] = barycentric2cartesian(u(i),v(i),w(i),vert);
    end
else
    y = eta1;
    z = eta2;
end


%%
for i = 1:size(N,2)
    figure
    surf(y,z,reshape(N(:,i),size(eta1,1),[]))
end
% pause
end
