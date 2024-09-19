function shp = shpBsplinesP3(n,x,~)
%% B-spline shape functions of order 3
% performs k-refinement automatically such that n is the total number of shape
% functions.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
p=3;

nAdd=n-(p+1);
xk=[-ones(1,p),linspace(-1,1,nAdd+2),ones(1,p)];


shp = element.shpBsplines(n,x,false,xk,p);


