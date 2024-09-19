function  [obj, msh] = createMesh(obj,msh,~)
%% cross-section mesh of full cylinder (y,z)

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
%% interface
R = obj.R;
n = obj.n;

%% vertices
theta = linspace(0,2*pi,n+1);                                               % angles of vertices
theta = theta(1:end-1);                                                     % last one is the same as first
theta = theta - pi/2 + pi/n;

coord = zeros(n+1,3);                                                       % coordinate table
coord(1:n,2)  = cos(theta)*R;                                               % horizontal coordinate y
coord(1:n,3)  = sin(theta)*R;                                               % vertical coordinate z

if n < 5
    connectivity = 1:n;
    coord(end,:) = [];
else
    connectivity = [ones(1,n)*(n+1);(1:n);(1:n)+1].';
    connectivity(end) = 1;
end

msh.connectivity = connectivity;
msh.coord  = coord;
msh.coordV = coord;


end
