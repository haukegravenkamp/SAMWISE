
function [u,v,w,J]=local2Barycentric(eta1,eta2)

u=(eta1+1)/2;                                                               % interpret eta1 as u
v=(eta2+1)/2;
v=v.*(1-u);

w=1-u-v;

udeta1=1/2;
udeta2=0;
vdeta1=-(eta2+1)/4;
vdeta2=1/2-(eta1+1)/4;
wdeta1=-udeta1-vdeta1;
wdeta2=-udeta2-vdeta2;

J=[udeta1, vdeta1 wdeta1; udeta2 vdeta2 wdeta2];                                         % Jacobian

end
%   2012-2024 Hauke Gravenkamp, gravenkamp.research@gmail.com
 
