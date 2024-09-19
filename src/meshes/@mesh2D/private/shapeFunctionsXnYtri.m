
function [shp,shpDeta1,shpDeta2] = shapeFunctionsXnYtri(p,eta12,nModes,nDom,intTab,iT,proriolOrders,iV)

%%
pDom = p(1:2);
pSid = p(3:end);

%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
indE = indicesCornersEdges(p,3);                                            % indices of edge modes
indD = (max(indE(:))+1):nModes;                                             % indices of domain modes

[u,v,w,J] = local2Barycentric(eta12(:,1),eta12(:,2));                       % get Barycentric coordinates


%% COMPUTE SHAPE FUNCTIONS ALONG THE SIDES

% functions f,g,h see Eq. 2.105 in "Precursors..."
[f_v,f_vD] = edgeShp(v  ,pSid(1),nModes,nonnan(indE(1,:)),intTab);          % f(v), and d/dv f(v)
[f_u,f_uD] = edgeShp(1-u,pSid(1),nModes,nonnan(indE(1,:)),intTab);          % f(1-u), and d/d(1-u) f(1-u)
[f_0,~]    = edgeShp(0  ,pSid(1),nModes,nonnan(indE(1,:)),intTab);          % f(0), and f'(0)

[g_w,g_wD] = edgeShp(w  ,pSid(2),nModes,nonnan(indE(2,:)),intTab);          % g(w), and d/dw g(w)
[g_v,g_vD] = edgeShp(1-v,pSid(2),nModes,nonnan(indE(2,:)),intTab);          % g(1-v), and d/d(1-v) g(1-v)
[g_0,~]   =  edgeShp(0  ,pSid(2),nModes,nonnan(indE(2,:)),intTab);          % g(0), and g'(0)

[h_u,h_uD] = edgeShp(u  ,pSid(3),nModes,nonnan(indE(3,:)),intTab);                                                % h(u), and d/du h(v)
[h_w,h_wD] = edgeShp(1-w,pSid(3),nModes,nonnan(indE(3,:)),intTab);                                              % h(1-w), and d/d(1-w) h(1-w)
[h_0,~]    = edgeShp(0  ,pSid(3),nModes,nonnan(indE(3,:)),intTab);                                                % h(0), and h'(0)


% [f_u,f_uD]=edgeShp(obj,1-u,1);                                              % f(1-u), and d/d(1-u) f(1-u)
% [f_0,~]   =edgeShp(obj,0,1);                                                % f(0), and f'(0)

% [g_w,g_wD]=edgeShp(obj,w,2);                                             % g(w), and d/dw g(w)
% [g_v,g_vD]=edgeShp(obj,1-v,2);                                              % g(1-v), and d/d(1-v) g(1-v)
% [g_0,~]   =edgeShp(obj,0,2);                                                % g(0), and g'(0)

% [h_u,h_uD]=edgeShp(obj,u,3);                                                % h(u), and d/du h(v)
% [h_w,h_wD]=edgeShp(obj,1-w,3);                                              % h(1-w), and d/d(1-w) h(1-w)
% [h_0,~]   =edgeShp(obj,0,3);                                                % h(0), and h'(0)


shp=1/2* (...
    v*g_w/(1-w) + ...
    u*h_w/(1-w) + ...
    w*h_u/(1-u) + ...
    v*f_u/(1-u) + ...
    u*f_v/(1-v) + ...
    w*g_v/(1-v) - ...
    u*f_0 - v*g_0 - w*h_0);

shpDu = 1/2* (                     ...
    + 0                            ...
    +   h_w /(1-w)                 ...
    + w*h_uD/(1-u) + w*h_u/(1-u)^2 ...
    - v*f_uD/(1-u) + v*f_u/(1-u)^2 ...
    +   f_v/(1-v)                  ...
    + 0                            ...
    - f_0 );

shpDv = 1/2* (                      ...
    +   g_w /(1-w)                  ...
    + 0                             ...
    + 0                             ...
    +   f_u /(1-u)                  ...
    + u*f_vD/(1-v) + u*f_v/(1-v)^2  ...
    - w*g_vD/(1-v) + w*g_v/(1-v)^2  ...
    - g_0);

shpDw = 1/2* (                     ...
    + v*g_wD/(1-w) + v*g_w/(1-w)^2 ...
    - u*h_wD/(1-w) + u*h_w/(1-w)^2 ...
    +   h_u /(1-u)                 ...
    + 0                            ...
    + 0                            ...
    +   g_v /(1-v)                 ...
    - h_0);


shpD=J*[shpDu;shpDv;shpDw];

shpDeta1=shpD(1,:);
shpDeta2=shpD(2,:);

[shp,shpDeta1,shpDeta2]=interiorModesTri(shp,shpDeta1,shpDeta2,eta12(:,1),eta12(:,2),indD,p(1),proriolOrders,iV);

   

%% Apply blending to create nodal shape functions
if ~isempty(iT)
    shp = shp*iT;
    shpDeta1 = shpDeta1*iT;
    shpDeta2 = shpDeta2*iT;
end

end
