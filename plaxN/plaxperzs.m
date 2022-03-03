%**********************************************************************

function [ipdaC] = plaxperzs(mF,mFB,GDt,ipda0,ipdaB,ipdaC);

%======================================================================
% Get integration point data from ipda0, ipdaB and ipdaC
% Variable names from 'old' implementation ('avpcpls') are mentioned.
% Following variables from 'avpcpls' are not used any more :
% F, F0, thk0, t33, t330

% Material parameters are taken from ipda0.

E    = ipda0(6); 
nu   = ipda0(7); 
thk0 = ipda0(3);                   % thick0;
Gty  = ipda0(10);                  % GsY;
h1   = ipda0(11);                  % hard;
Gg   = ipda0(8);  
N    = ipda0(9);
Gm   = E/(2*(1+nu));               % mu;
Gl   = (E*nu)/((1+nu)*(1-2*nu));   % lambda;
h2   = ipda0(12);                  % AAA;
h3   = ipda0(13);                  % BBB;
h4   = ipda0(14);                  % CCC;
h5   = ipda0(15);                  % DDD;

% Other data (begin increment and current) are taken from database.

ccGs  = ipdaC(60:64)';             % Gs;
ccGt  = ipdaC(65:69)';             % Gt';
ccGtB = ipdaB(65:69)';             % Gt0';
GDGl  = ipdaC(76);                 % Glinc;
Gk    = ipdaC(77);        
GkB   = ipdaB(77);                 % Gk0;
GDGk  = ipdaC(79);                 % Gkinc;
DGk   = ipdaC(80);                 % Gkd;
Gz    = ipdaC(81);                 % Y;

%======================================================================
% Initialize some variables
% mI   = unity matrix
% mmI  = 4th-order unity marix
% mmIs = symmetric 4th-order unity matrix

mI   = eye(3);     
mmI  = m2mm(mI,5);
mmI  = mmI(:,[1 2 3 5 4]);
mmIs = 1/2*(mmI + mmI(:,[1 2 3 5 4]));
ccI  = m2cc(mI,5);

mGs = zeros(3);
mGs(1,1) = ccGs(1); mGs(2,2) = ccGs(2); mGs(3,3) = ccGs(3);
mGs(1,2) = ccGs(4); mGs(2,1) = ccGs(5);

mGt = zeros(3);
mGt(1,1) = ccGt(1); mGt(2,2) = ccGt(2); mGt(3,3) = ccGt(3);
mGt(1,2) = ccGt(4); mGt(2,1) = ccGt(5);

mGtB = zeros(3);
mGtB(1,1) = ccGtB(1); mGtB(2,2) = ccGtB(2); mGtB(3,3) = ccGtB(3);
mGtB(1,2) = ccGtB(4); mGtB(2,1) = ccGtB(5);

J    = det(mF); 
mFi  = inv(mF);
mFBi = inv(mFB);
mFn  = mF * mFBi;
mFni = inv(mFn);
men  = 1/2*(mI - mFni' * mFni);
ccen = m2cc(men,5);

%======================================================================
% elastic material matrix

mmH  = 2*(Gm - Gl*log(J)) .* mmIs + Gl*ccI*ccI';

%----------------------------------------------------------------------
% Incremental deformation is assumed to be elastic and the elastic
% predictor stress (trial stress) is calculated.
% The flow criterion is evaluated.

mA   = mFn * mGtB * mFn';
ccA  = m2cc(mA,5);

ccGte  = ccA + mmH(:,[1 2 3 5 4]) * ccen;
ccGted = ccGte - 1/3 * sum(ccGte(1:3)) * ccI;
Gteq   = sqrt( 3/2 * ccGted' * ccGted([1 2 3 5 4]) );
cca    = 3/2 * 1/Gteq * ccGted;

Gz = Gty + h1*Gk + h2*Gk^2 + h3*Gk^3 + h4*Gk^4 + h5*Gk^7;
F  = Gteq - Gz;

%----------------------------------------------------------------------
if (F>=0)
%----------------------------------------------------------------------
% viscoplastic deformation

% initial update of stress

Gf    = (F/Gty)^N;
Gf    = 1/2 * (Gf + abs(Gf));
GDGl  = GDt * Gg * Gf;

ccGt  = ccGte - GDGl * mmH(:,[1 2 3 5 4]) * cca;
ccGtd = ccGt - 1/3 * sum(ccGt(1:3)) * ccI;
Gteq  = sqrt( 3/2 * ccGtd' * ccGtd([1 2 3 5 4]) );
cca   = 3/2 * 1/Gteq * ccGtd;

DGk   = Gg * Gf * sqrt( 2/3 * cca' * cca );
GDGk  = GDt * DGk;
Gk    = GkB + GDGk;

Gz    = Gty + h1*Gk + h2*Gk^2 + h3*Gk^3 + h4*Gk^4 + h5*Gk^7;
F     = Gteq - Gz;
Gf    = (F/Gty)^N;
Gf    = 1/2 * (Gf + abs(Gf));
GfF   = N * (F/Gty)^(N-1) * 1/Gty;
GfF   = 1/2 * (GfF + abs(GfF));
FGk   = - h1 - 2*h2*Gk - 3*h3*Gk^2 - 4*h4*Gk^3 - 7*h5*Gk^6;
mmb   = 1/Gteq * ( - cca * cca' + 3/2 * mmIs - 1/2 * ccI * ccI');

% initial system

mmR   = mmI + GDGl * mmH(:,[1 2 3 5 4]) * mmb;
cct   = mmH * cca([1 2 3 5 4]);
ccu   = - GDt * Gg * GfF * cca;
v     = 1 - GDt * Gg * GfF * FGk;

ccs1  = ccGt - ccGte + GDGl * mmH * cca([1 2 3 5 4]);
s2    = GDGl - GDt * Gg * Gf;

% start iteration loop

ii    = 0;
eps1  = 1e-8;
eps2  = 1e-8;
mii   = 1000;

while (ii==0 | nrm1>eps1 | nrm2>eps2)
%----------------------------------------------------------------------
ii = ii + 1;

if ii>=mii
  fprintf(1,'no convergence in stress-update algorithm\n');
  [norm(ccs1,2) s2 ii]
  break
end

% There are two ways to solve the iterative unknowns.
% In the first method *one* matrix is build and the unknowns are
% solved simultaneously.
% In the second method the unknowns are solved subsequently,
% as was also implemented in 'avpcpls'.
% Results of the two methods are the same.

% simultaneous solution

%LHM    = [mmR cct ; ccu' v];
%RHC    = [-ccs1 ; -s2];
%LHMi   = inv(LHM);
%sol    = LHMi * RHC;
%ccGdGt = sol(1:5);
%GdGl   = sol(6);

% subsequent solution

mmRi   = inv(mmR);
mmRi   = mmRi(:,[1 2 3 5 4]);
GdGl   = (-s2 + ccu'*mmRi*ccs1)/(v - ccu'*mmRi*cct);
ccGdGt = - mmRi * ( ccs1 + cct*GdGl );

% update stress and GDGl

ccGt  = ccGt + ccGdGt;
GDGl  = GDGl + GdGl;

% make system again

ccGtd = ccGt - 1/3 * sum(ccGt(1:3)) * ccI;
Gteq  = sqrt( 3/2 * ccGtd' * ccGtd([1 2 3 5 4]) );
cca   = 3/2 * 1/Gteq * ccGtd;

DGk   = Gg * Gf * sqrt( 2/3 * cca' * cca );
GDGk  = GDt * DGk;
Gk    = GkB + GDGk;

Gz    = Gty + h1*Gk + h2*Gk^2 + h3*Gk^3 + h4*Gk^4 + h5*Gk^7;
F     = Gteq - Gz;
Gf    = (F/Gty)^N;
Gf    = 1/2 * (Gf + abs(Gf));
GfF   = N * (F/Gty)^(N-1) * 1/Gty;
GfF   = 1/2 * (GfF + abs(GfF));
FGk   = - h1 - 2*h2*Gk - 3*h3*Gk^2 - 4*h4*Gk^3 - 7*h5*Gk^6;
mmb   = 1/Gteq * ( - cca * cca' + 3/2 * mmIs - 1/2 * ccI * ccI');

% new system

mmR   = mmI + GDGl * mmH(:,[1 2 3 5 4]) * mmb;
cct   = mmH * cca([1 2 3 5 4]);
ccu   = - GDt * Gg * GfF * cca;
v     = 1 - GDt * Gg * GfF * FGk;

ccs1  = ccGt - ccGte + GDGl * mmH * cca([1 2 3 5 4]);
s2    = GDGl - GDt * Gg * Gf;

nrm1  = norm(ccs1,2);
nrm2  = abs(s2);

%----------------------------------------------------------------------
end;

%----------------------------------------------------------------------
else
%----------------------------------------------------------------------
% elastic deformation

 mGt   = Gm * ( mF * mF' - mI ) + Gl * log(J) * mI;
% ccGt  = m2cc(mGt,5);
 ccGt  = ccGte;

DGk   = 0;
ccD   = zeros(5,1);
GDGk  = 0;
Gk    = GkB;

%----------------------------------------------------------------------
end;
%----------------------------------------------------------------------

ccGs = 1/J * ccGt;

%======================================================================
% update ipdaC

ipdaC(60:64) = ccGs;
ipdaC(65:69) = ccGt;
ipdaC(76)    = GDGl;
ipdaC(77)    = Gk;
ipdaC(78)    = ccGt(3);
ipdaC(79)    = GDGk;
ipdaC(80)    = DGk;
ipdaC(81)    = Gz;
%ipdaC(82)    = mis;

%**********************************************************************
