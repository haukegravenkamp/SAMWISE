function [x,y,J]=coordinateTransformTri(obj,eta1,eta2,~)
%% Transformation from local to Cartesian coordinates
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 

[u,v,~,J1]=local2Barycentric(eta1,eta2);

if obj.curved                                                               % element has curved edges
    % TO DO
else
    
    
    % vertex coordinates
    xn=obj.vCoord(obj.masterNodes,1);
    yn=obj.vCoord(obj.masterNodes,2);
    
    [x,y,J2] = barycentric2cartesian(u,v,1-u-v,[xn,yn]);

    
    
end
J=J1(1:2,1:2)*J2;

end
