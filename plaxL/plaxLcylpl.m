%======================================================================
% Thick-walled cylinder with internal pressure; plane stress;
%----------------------------------------------------------------------
delete('loadincr.m','savedata.m'); 

a = 0.25; b = 0.5; h = 0.5; pin = 100e6; pex = 0; E = 250e9; Gn = 0.33; 

elda  = [ 1 11 0.5 0 0 E Gn 0 0 0 0 0 0 0 0 ];

pp = [ 
 17 1 0; 34 1 0; 51 1 0; 68 1 0; 85 1 0; 102 1 0; 119 1 0; 136 1 0; 153 1 0;
  1 2 0; 18 2 0; 35 2 0; 52 2 0; 69 2 0;  86 2 0; 103 2 0; 120 2 0; 137 2 0 ];

fp = [ 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 ]; lfp = length(fp);
for i=1:lfp
  np = fp(i); xxx = crd0(np,1); yyy = crd0(np,2); 
  rrr = sqrt(xxx^2+yyy^2); cos = xxx/rrr; sin = yyy/rrr; 
  pf(2*(i-1)+1,:) = [ np 1 pin*h*pi/2*a/(lfp-1) * cos ];
  pf(2*(i-1)+2,:) = [ np 2 pin*h*pi/2*a/(lfp-1) * sin ];
end;

pf(1:2,3) = pf(1:2,3)./2; pf((2*lfp-1):2*lfp,3) = pf((2*lfp-1):2*lfp,3)/2;

%======================================================================