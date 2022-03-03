%**********************************************************************

function [ipdaC,ipsmC] = plaxviels(GDt,ipda0,ipdaB,ipdaC,ipsmB,eldalivi);

einf = ipda0(6);
nuu  = ipda0(7);
nmo  = ipda0(8);
Kinf = einf/(3*(1-2*nuu));
Ginf = einf/(2*(1+nuu));

ccGe = [ipdaC(17:20) ipdaC(20)];
ccGDGe = ccGe - [ipdaB(17:20) ipdaB(20)];

fac = 1/3;
mmAh = [ fac fac fac 0 0 
         fac fac fac 0 0 
         fac fac fac 0 0 
           0   0   0 0 0 
           0   0   0 0 0 ];
mmAd = [ 2*fac  -fac  -fac 0 0 
          -fac 2*fac  -fac 0 0 
          -fac  -fac 2*fac 0 0 
             0     0     0 1 0
             0     0     0 0 1 ];

CC = (3*Kinf*mmAh + 2*Ginf*mmAd)*ccGe';

SS = zeros(5,1);

for i=1:nmo
%----------------------------------------------------------------------
Ei  = eldalivi(i,1);
Gti = eldalivi(i,2);
Ki  = Ei/(3*(1-2*nuu));
Gi  = Ei/(2*(1+nuu));

ee = exp(-GDt/Gti);
ccGsiB = ipsmB(5*(i-1)+1:5*i);
pp = (Gti/GDt)*(1-ee);
ccGsi = ee*ccGsiB + (pp*(3*Ki*mmAh + 2*Gi*mmAd)*ccGDGe')';
ipsmC(5*(i-1)+1:5*i) = ccGsi;
SS = SS + ccGsi';
%----------------------------------------------------------------------
end;

ccGs = (CC + SS)';
ipdaC(21:24) = ccGs(1:4);

%**********************************************************************
