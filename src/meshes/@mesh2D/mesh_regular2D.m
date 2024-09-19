function [v,e]=mesh_regular2D(n1,n2,L1,L2)

if nargin<1
    n1=5;                                   % number of elements in first directions
end
if nargin<2
    n2=4;                                   % number of elements in second direction
end
if nargin<3
    L1=1;
end
if nargin<4
    L2=1;
end

v1=linspace(0,L1,n1+1);
v2=linspace(0,L2,n2+1);

[V1,V2]=meshgrid(v1,v2);
v=[V1(:),V2(:)];

i=(1:n1*(n2+1))';
e=[i+1,i,i+n2+1,i+n2+2];
e(n2+1:n2+1:end,:)=[];
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
