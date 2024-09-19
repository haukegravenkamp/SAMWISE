function shp = shpHier(nNodes,x,~)
%% SHPHIER
% Hierarchical shape functions
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
if nNodes<12                                                                % use analytical expressions for efficiency
    
    N0=1/2*(1-x);                                                           % linear shape functions
    N1=1/2*(1+x);
    
    shp(1,1)=N0;
    shp(1,2)=N1;
    shp(1,3)=(-2*sqrt(3/2));
    shp(1,4)=(-2*sqrt(5/2)*x);
    shp(1,5)=(-1/2*sqrt(7/2)*(5*x.^2-1));
    shp(1,6)=(-1/2*sqrt(9/2)*(7*x.^2-3).*x);
    shp(1,7)=(-1/4*sqrt(11/2)*(21*x.^4-14*x.^2+1));
    shp(1,8)=(-1/4*sqrt(13/2)*(33*x.^4-30*x.^2+5).*x);
    shp(1,9)=(-1/32*sqrt(15/2)*(429*x.^6-495*x.^4+135*x.^2 - 5));
    shp(1,10)=(-1/32*sqrt(17/2)*(715*x.^6-1001*x.^4+385*x.^2 - 35).*x);
    shp(1,11)=(-1/64*sqrt(19/2)*(2431*x.^8-4004*x.^6+2002*x.^4 - 308.*x.^2+7));
    
    
    C=-1/2*1/2*(1+x)+1/2*(1-x)*1/2;
    shp(2,1)=-1/2;
    shp(2,2)=1/2;
    shp(2,3)=C*shp(1,3);
    shp(2,4)=C*shp(1,4) + N0.*N1.* (-2*sqrt(5/2));
    shp(2,5)=C*shp(1,5) + N0.*N1.* (-1/2*sqrt(7/2)*(2*5*x));
    shp(2,6)=C*shp(1,6) + N0.*N1.* (-1/2*sqrt(9/2)*(3*7*x.^2-3));
    shp(2,7)=C*shp(1,7) + N0.*N1.* (-1/4*sqrt(11/2)*(4*21*x.^3-2*14*x));
    shp(2,8)=C*shp(1,8) + N0.*N1.* (-1/4*sqrt(13/2)*(5*33*x.^4-3*30*x.^2+5));
    shp(2,9)=C*shp(1,9) + N0.*N1.* (-1/32*sqrt(15/2)*(6*429*x.^5-4*495*x.^3 + 2*135*x ));
    shp(2,10)=C*shp(1,10) + N0.*N1.* (-1/32*sqrt(17/2)*(7*715*x.^6-5*1001*x.^4 + 3*385*x.^2 - 35));
    shp(2,11)=C*shp(1,11) + N0.*N1.* (-1/64*sqrt(19/2)*(8*2431*x.^7-6*4004*x.^5+4*2002*x.^3 - 2*308.*x));
    
    shp(1,3:end)=N0.*N1.*shp(1,3:end);
    
    shp=shp(:,1:nNodes);
    
    
    shp=[shp(:,1),shp(:,3:end),shp(:,2)];
    
else
    
    N0=1/2*(1-x);
    N1=1/2*(1+x);
    
    ks=1:nNodes-2;                                                          % required orders of Legendre polynomials to integrate
    
    ki=0:nNodes-1;                                                          % orders of all Legendre polynomials involved
%     LP=legendreP(ki,x);                                                   % all Legendre polynomials
    LP=element.jacobiP(0,0,ki,x);                                           % all Jacobi polynomials
    LPk=LP(2:end-1);                                                        % L(k)
    LPkp1=LP(3:end);                                                        % L(k+1)
    LPkm1=LP(1:end-2);                                                      % L(k-1)
    
    Lob=(LPkp1-LPkm1)./(2*ks+1);
    Lob=sqrt((2*ks+1)/2).*Lob;                                              % normalize
    
    LobD = sqrt((2*ks+1)/2) .* LPk;
    
    shp=[N0 Lob N1; -0.5 LobD 0.5];
    
end




