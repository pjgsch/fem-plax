%**********************************************************************

function [em,ef,Gdn,Gdt,eidaC] = plaxczelem1(em,ef,eidaC,d0,Gfn,Gft,Tnmax,Ttmax,ec0,ec,ip,gip,Tpe);

% em    = element stiffness matrix
% ef    = element internal force column
% Gdn   = delta_n : characteristic normal opening
% Gdt   = delta_t : characteristic tangential opening
% eidaC = data base array with current values
%
% d0    = initial 'depth' of cz element
% Gfn   = fie_n : work-of-separation in normal direction
% Gft   = fie_t : work-of-separation in tangential direction
% Tnmax = T_{n,max} : maximum normal force
% Ttmax = T_{t,max} : maximum tangential force
% ec0   = initial element nodal coordinates 
% ec    = current element nodal coordinates
% ip    = integration point number
% gip   = global integration point number 
% Tpe   = total displacement of element nodes

Gdn = Gfn/(Tnmax*exp(1));
Gdt = Gft/(Ttmax*sqrt(0.5*exp(1)));

x10 = ec0(1,1); x20 = ec0(2,1); x30 = ec0(3,1); x40 = ec0(4,1);
y10 = ec0(1,2); y20 = ec0(2,2); y30 = ec0(3,2); y40 = ec0(4,2);
xA0 = (x10+x40)/2; xB0 = (x20+x30)/2;
yA0 = (y10+y40)/2; yB0 = (y20+y30)/2;
l0  = sqrt( (yB0-yA0)^2 + (xB0-xA0)^2 );

x1 = ec(1,1); x2 = ec(2,1); x3 = ec(3,1); x4 = ec(4,1);
y1 = ec(1,2); y2 = ec(2,2); y3 = ec(3,2); y4 = ec(4,2);
xA = (x1+x4)/2; xB = (x2+x3)/2; yA = (y1+y4)/2; yB = (y2+y3)/2; 
%xA = x4; xB = x3; yA = y4; yB = y3;   % plane AB on bottom
%xA = x1; xB = x2; yA = y1; yB = y2;   % plane AB on top
l  = sqrt( (yB-yA)^2 + (xB-xA)^2 ); 
c  = (xB-xA)/l; s = (yB-yA)/l;

R = [  c  s  0  0  0  0  0  0  
      -s  c  0  0  0  0  0  0 
       0  0  c  s  0  0  0  0 
       0  0 -s  c  0  0  0  0 
       0  0  0  0  c  s  0  0 
       0  0  0  0 -s  c  0  0 
       0  0  0  0  0  0  c  s 
       0  0  0  0  0  0 -s  c ];

P = [ -1  0  0  0  0  0  1  0 
       0 -1  0  0  0  0  0  1 
       0  0 -1  0  1  0  0  0 
       0  0  0 -1  0  1  0  0 ];

GDuAB = P * R * Tpe;

if     ip==1, ksi = -1/sqrt(3); ifac = 1;
elseif ip==2, ksi =  1/sqrt(3); ifac = 1;
end;

NN = 0.5 * [ 1-ksi 0 1+ksi 0 ; 0 1-ksi 0 1+ksi ];
GDu = NN * GDuAB;
GDut = GDu(1); GDun = GDu(2);

an = Gfn/Gdn;  at = Gft/Gdt;
un = GDun/Gdn; ut = GDut/Gdt;

TTn = an * un * exp(-un) * exp(-ut^2);
TTt = 2*at * ut * (1+un) * exp(-un) * exp(-ut^2);
TT  = [ TTt  TTn ]';
MM(1,1) = ( 2*at * (1+un)  * (1-2*ut^2) * exp(-un) * exp(-ut^2))/Gdt; % DtDt
MM(1,2) = (-2*at * un * ut * exp(-un) * exp(-ut^2)             )/Gdn; % DtDn
MM(2,1) = (-2*an * un * ut * exp(-un) * exp(-ut^2)             )/Gdt; % DnDt
MM(2,2) = (   an * (1-un)  * exp(-un) * exp(-ut^2)             )/Gdn; % DnDn
        
eidaC(gip,20) = GDun;
eidaC(gip,21) = GDut;
eidaC(gip,22) = TTn;
eidaC(gip,23) = TTt;

em = em + (d0*l0/2) * R' * P' * NN' * MM * NN * ifac * P * R ;
ef = ef + (d0*l0/2) * R' * P' * NN' * TT * ifac ;

%**********************************************************************
