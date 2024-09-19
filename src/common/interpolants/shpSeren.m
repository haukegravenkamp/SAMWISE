function shp = shpSeren(nNodes,xi,~)
%% Serendipity shape functions along the edges
%
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
shp=zeros(2,nNodes);                                                        % allocate shape functions
p=nNodes-1;                                                                 % element order
eta=-1;                                                                     % use the 2D shape functions at eta=-1

eta2 = eta.*eta;
eta3 = eta.*eta.*eta;
eta4 = eta.*eta.*eta.*eta;
etaxi = eta.*xi;
eta2xi = eta.*eta.*xi;
etaxi2 = eta.*xi.*xi;
xi2 = xi.*xi;
xi3 = xi.*xi.*xi;
xi4 = xi.*xi.*xi.*xi;

switch p                                                                    % switch element orders
    case 1
        shp(1,1) = 1/4*(1-xi).*(1-eta);
        shp(1,2) = 1/4*(1+xi).*(1-eta);
        shp(2,1) = -1/4*(1-eta);
        shp(2,2) = 1/4*(1-eta);
    case 2
        shp(1,1) = -(1-xi).*(1-eta).*(1+xi+eta)/4;
        shp(1,3) = -(1+xi).*(1-eta).*(1-xi+eta)/4;
        shp(1,2) = (1-xi.*xi).*(1-eta)/2;
        
        shp(2,1) = -(-(1-eta).*(1+xi+eta)+(1-xi).*(1-eta))/4;
        shp(2,3) = -( (1-eta).*(1-xi+eta)-(1+xi).*(1-eta))/4;
        shp(2,2) = -xi.*(1-eta);
    case 3
        shp(1,1) =  1/32*(eta-1).*(xi-1).*(9*eta.*eta+9*xi.*xi-10);
        shp(1,4) = -1/32*(eta-1).*(xi+1).*(9*eta.*eta+9*xi.*xi-10);
        shp(1,2) =  9/32*(eta-1).*(-3*xi.*xi.*xi+xi.*xi+3*xi-1);
        shp(1,3) = -9/32*(eta-1).*(-3*xi.*xi.*xi-xi.*xi+3*xi+1);
        
        
        shp(2,1) = -1/32*(eta-1).*(-9*eta.*eta-27*xi.*xi+18*xi+10);
        shp(2,4) = -1/32*(eta-1).*(9*eta.*eta+27*xi.*xi+18*xi-10);
        shp(2,2) =  9/32*(eta-1).*(-9*xi.*xi+2*xi+3);
        shp(2,3) =  9/32*(eta-1).*(9*xi.*xi+2*xi-3);
        
    case 4
        shp(1,1) =  1/12*(eta-1).*(xi-1).*(-4*eta3+3*etaxi+4*eta-4*xi3+4*xi);
        shp(1,5) =  1/12*(eta-1).*(xi+1).*(4*eta3+3*etaxi-4*eta-4*xi3+4*xi);
        shp(1,2) = -2/3*xi.*(eta-1).*(-2*xi3+xi2+2*xi-1);
        shp(1,3) = -1/2*(xi2-1).*(eta-1).*(4*xi2+eta);
        shp(1,4) = -2/3*xi.*(eta-1).*(-2*xi3-xi2+2*xi+1);
        
        shp(2,1) =  1/12*(eta-1).*(-4*eta3+6*etaxi+eta-16*xi3+12*xi2+8*xi-4);
        shp(2,5) =  1/12*(eta-1).*(4*eta3+6*etaxi-eta-16*xi3-12*xi2+8*xi+4);
        shp(2,2) = -2/3*(eta-1).*(-8*xi3+3*xi2+4*xi-1);
        shp(2,3) = -xi.*(eta-1).*(8*xi2+eta-4);
        shp(2,4) = -2/3*(eta-1).*(-8*xi3-3*xi2+4*xi+1);
    case 5
        shp(1,1) = -1/24576*(eta-1).*(xi-1).*(-10000*eta4+13462*eta2xi+17462*eta2+17616*etaxi2+20193*etaxi+2577*eta-10000*xi4+21616*xi2+6731*xi-5029);
        shp(1,6) = -1/24576*(eta-1).*(xi+1).*(10000*eta4+13462*eta2xi-17462*eta2-17616*etaxi2+20193*etaxi-2577*eta+10000*xi4-21616*xi2+6731*xi+5029);
        shp(1,2) =  25/6144*(xi2-1).*(eta-1).*(134*eta2+210*etaxi+75*eta-500*xi3+300*xi2+230*xi-71);
        shp(1,3) = -25/3072*(xi2-1).*(eta-1).*(2*eta2-110*etaxi+25*eta-500*xi3+100*xi2+70*xi-13);
        shp(1,4) = -25/3072*(xi2-1).*(eta-1).*(2*eta2+110*etaxi+25*eta+500*xi3+100*xi2-70*xi-13);
        shp(1,5) =  25/6144*(xi2-1).*(eta-1).*(134*eta2-210*etaxi+75*eta+500*xi3+300*xi2-230*xi-71);
        
        shp(2,1) = -1/12288*(eta-1).*(-5000*eta4+13462*eta2xi+2000*eta2+26424*etaxi2+2577*etaxi-8808*eta-25000*xi4+20000*xi3+32424*xi2-14885*xi-5880);
        shp(2,6) = -1/12288*(eta-1).*(5000*eta4+13462*eta2xi-2000*eta2-26424*etaxi2+2577*etaxi+8808*eta+25000*xi4+20000*xi3-32424*xi2-14885*xi+5880);
        shp(2,2) =  25/3072*(eta-1).*(134*eta2xi+315*etaxi2+75*etaxi-105*eta-1250*xi4+600*xi3+1095*xi2-371*xi-115);
        shp(2,3) = -25/1536*(eta-1).*(2*eta2xi-165*etaxi2+25*etaxi+55*eta-1250*xi4+200*xi3+855*xi2-113*xi-35);
        shp(2,4) = -25/1536*(eta-1).*(2*eta2xi+165*etaxi2+25*etaxi-55*eta+1250*xi4+200*xi3-855*xi2-113*xi+35);
        shp(2,5) =  25/3072*(eta-1).*(134*eta2xi-315*etaxi2+75*etaxi+105*eta+1250*xi4+600*xi3-1095*xi2-371*xi+115);
end






