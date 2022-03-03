%**********************************************************************

function [mmM,mmGS] = plaxvielm(GDt,ipda0,eldalivi,it);

einf = ipda0(6);
nuu  = ipda0(7);
nmo  = ipda0(8);
Kinf = einf/(3*(1-2*nuu));
Ginf = einf/(2*(1+nuu));

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

mmM = 3*Kinf*mmAh + 2*Ginf*mmAd;

for i=1:nmo
%----------------------------------------------------------------------
Ei  = eldalivi(i,1);
Gti = eldalivi(i,2);
Ki  = Ei/(3*(1-2*nuu));
Gi  = Ei/(2*(1+nuu));

ee = exp(-GDt/Gti);
pp = (Gti/GDt)*(1-ee);
if it==0, pp=1; end;

mmM = mmM + pp*(3*Ki*mmAh + 2*Gi*mmAd);
%----------------------------------------------------------------------
end;

mmGS = zeros(5,5);

%**********************************************************************
