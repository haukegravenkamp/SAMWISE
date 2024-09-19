function [shp,wn] = shpMLS(nNodes,eta,~)
%% Moving least squares shape functions
% He, Y., Yang, H., & Deeks, A. J. (2012).
% An element-free Galerkin (EFG) scaled boundary method.
% Finite Elements in Analysis and Design, 62, 28–36.
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
type=2;                                                                     % type of weight function

if type==0                                                                  % constant weight function (leads to standard fem)
    order=nNodes-1;                                                         % order of polynomial basis
else                                                                        % spline weight function
    if nNodes<4
        order=1;
    else
        order=2;                                                                % order of polynomial basis
    end
end


beta=nNodes-1;                                                                     % parameter controlling the support size
% beta=3;

s=(eta+1)/2;                                                                % element parametrisation (0->1)

ds=1/(nNodes-1);                                                            % node distance
ri=beta*ds;                                                                 % support size

sn=linspace(0,1,nNodes);                                                    % nodal positions

[phi,dphi,wn]=computePhi(s,sn);

T=zeros(nNodes);

for iN=1:nNodes
    T(iN,:)=computePhi(sn(iN),sn);
end
iT=T\eye(nNodes);
phi=phi*iT;
dphi=dphi*iT;


dphi=1/2*dphi;                                                              % coordinate change s -> eta
shp=[phi;dphi];





%% COMPUTE BASIS FUNCTIONS
    function [phi,dphi,wn]=computePhi(s,sn)
        
        
        
        [wn,dwn]=computeWn(s,sn,type);                                                   % weight function and derivative
        
        
        % compute matrices A,B according to Eqs. (31-33)
        A=zeros(order+1);                                                           % initialize matrix A
        B=zeros(order+1,nNodes);                                                    % initialize matrix B
        dA=zeros(order+1);                                                          % initialize derivative of matrix A
        dB=zeros(order+1,nNodes);                                                   % initialize derivative of matrix B
        
        for i=1:nNodes                                                              % loop over nodes
            
            ps=computePs(sn(i),order);                                              % monomial basis
            A =  A +  wn(i)*(ps*ps');                                               % contribution to matrix A
            dA= dA + dwn(i)*(ps*ps');                                               % contribution to derivative
            B(:,i)  =  wn(i)*ps;                                                    % contribution to matrix B
            dB(:,i) = dwn(i)*ps;                                                    % contribution to derivative
            
        end
        
        [ps,dps]=computePs(s,order);                                                % monomial basis and derivative
        
        iA=A\eye(size(A));                                                          % inverse of A
        phi=ps'*(iA*B);
        
        diA=-iA*dA*iA;
        dphi=dps'*(iA*B) + ps'*(diA*B + iA*dB);
        
    end

%% SPLINE WEIGHT FUNCTION
    function [wn,dwn]=computeWn(s,sn,type)
        d =abs(s-sn);                                                       % distance from support
        r=abs(s-sn)/ri;                                                     % relative distance from support
        signD=sign(s-sn);                                                   % sign of d
        dR=1/ri*signD;                                                      % derivative of r wrt s
        dD=signD;                                                           % derivative of d wrt s
        switch type                                                         % types of weight functions
            case 0                                                          % constant
                wn  = ones(size(sn));
                dwn = zeros(size(sn));
            case 1                                                          % splines as used by He
                
                wn  =     (1-   6*r.^2 +   8*r.^3 -   3*r.^4).*(r<=1);      % weight function
                dwn = dR.*(0- 2*6*r.^1 + 3*8*r.^2 - 4*3*r.^3).*(r<=1);      % weight function derivative
                
            case 2                                                          % exponential weight function used by Zhu et al.
                ri=1;
                ci=ri/4;                                                    % value recommended by Zhu
                
                k=1;
                eR=exp(-(ri/ci)^(2*k));
                eD=exp(-(d/ci).^(2*k));
                wn=(eD - eR)/(1 - eR);
                deD = eD .*  (-2*k*(d/ci).^(2*k-1))/ci.*dD;
                dwn = deD/(1 - eR);
        end
        
    end

%% MONOMIAL BASIS
    function [ps,dps]=computePs(s,order)
        
        expo=(0:order)';                                                    % exponents
        ps=s.^expo;                                                         % evaluate basis
        
        dexpo=expo-1;                                                       % exponents of derivative
        dps=expo.*(s.^dexpo);                                               % evaluate derivative
        dps(1)=0;                                                           % remove first entry (constant term)
        
    end

end









