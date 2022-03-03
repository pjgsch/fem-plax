%**********************************************************************
%  This is a rigid rotation of a square
%  All nodal displacements are prescribed
%----------------------------------------------------------------------
clear all;delete('loadincr.m');delete('updaelda.m');delete('savedata.m');

neip  = 4;

lok   = [ 3 1 1 2 3 4 ];
crd0  = [ 0 0; 1 0; 1 1; 0 1 ];
elda  = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];

phi0 = 10/180*pi;
pp = [ 1 1 0; 1 2 0; ];
pp = [ pp; 2 1 -(1-cos(phi0)); 2 2 sin(phi0); ...
           4 1 -sin(phi0); 4 2 -(1-cos(phi0)); ...
           3 1 sqrt(2)*(cos(pi/4+phi0))-1; ...
           3 2 sqrt(2)*(sin(pi/4+phi0))-1];

nic    = 72;
nl     = 1;
mit    = 50;
deltat = 5;
cnm    = 5;
ccr    = 0.1;
vrs    = 1;
fat    = 0;

loin = fopen('loadincr.m','w');
fprintf(loin,'if ic>=1,\n');
fprintf(loin,'phi=ic*phi0;\n');
fprintf(loin,'x2 = cos(phi); y2 = sin(phi); \n');
fprintf(loin,'x3 = sqrt(2)*(cos(pi/4+phi)); y3 = sqrt(2)*(sin(pi/4+phi)); \n');
fprintf(loin,'x4 = -sin(phi); y4 = cos(phi); \n');
fprintf(loin,'u2 = x2-crd(2,1); v2 = y2-crd(2,2); \n');
fprintf(loin,'u3 = x3-crd(3,1); v3 = y3-crd(3,2); \n');
fprintf(loin,'u4 = x4-crd(4,1); v4 = y4-crd(4,2); \n');
fprintf(loin,'peC = [0; 0; u2; v2; u3; v3; u4; v4]; \n');
fprintf(loin,'end; \n');
fprintf(loin,'feC = ic*fe0; \n');
fclose(loin);

sada = fopen('savedata.m','w');
s1  = 'crd';
fprintf(sada,'save(savefile,s1); \n');
fprintf(sada,'Sic(ic) = ic; \n');
fprintf(sada,'Sfi(ic) = phi; \n');
fprintf(sada,'Spp(ic) = MDp(3,2);  \n');
fprintf(sada,'Sff(ic) = sqrt(Mfi(3,2)^2 + Mfi(3,1)^2); \n');
fprintf(sada,'Scx(ic) = crd(3,1); \n');
fclose(sada);

%**********************************************************************
