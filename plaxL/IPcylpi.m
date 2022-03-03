%**********************************************************************
%  axisymmetric : cylinder with internal pressure
%
%  nel = number of elements in radial direction
%  RRR = inner radius
%  TTT = wall thickness
%  LLL = height (length)
%  EMO = Young's modulus
%  POI = Poisson's ratio
%  PPP = internal pressure
%======================================================================
clear all;delete('loadincr.m');delete('updaelda.m');delete('savedata.m');
%======================================================================

% dimensions in meter

nel = 100; 
RRR = 0.25; TTT = 0.25; LLL = 0.5; 
EMO = 250e9; POI = 0.33;
PPP = 100e6;

neip  = 4;

base = [-1 1 2 0];
for i = 1:nel, lok(i,1:4) = base + i*2*[1 1 1 1]; end;
lok = [10.*ones(nel,1) ones(nel,1) lok ];

crd0 = [RRR 0; RRR LLL];
dx = TTT/nel;
for i=1:nel
  x = RRR + i*dx; yb = 0.0; yt = LLL;
  crd0 = [crd0; [x yb]; [x yt]];
end;

elda  = [10 1 1 0 0 EMO POI 0 0 0 0 0 0 0 0];

for i=1:nel+1, botn(i) = 2*(i-1)+1; topn(i) = 2*i; end;
pp    = [botn' 2*ones(nel+1,1) zeros(nel+1,1)];
%pp = [pp; topn' 2*ones(nel+1,1) zeros(nel+1,1)];

pl =  [topn(2:nel+1)' 2*ones(nel,1)];
pr =  [2 2];
lim = ones(nel,1);

FFF = PPP*2*pi*RRR*LLL;
pf  = [1 1 FFF/2; 2 1 FFF/2 ];

%vrs    = 3;

loin = fopen('loadincr.m','w');
fprintf(loin,'peC = pe0; \n');
fprintf(loin,'feC = fe0; \n');
fclose(loin);

sada = fopen('savedata.m','w');
fprintf(sada,'Sic = ic; \n');
fclose(sada);


%**********************************************************************
