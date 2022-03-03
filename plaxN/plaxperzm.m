%**********************************************************************

function [mmM,mmGS] = plaxperzm(mF,mFB,GDt,ipda0,ipdaB,ipdaC);

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

ccFBi    = m2cc(mFBi,5);
ccFn     = m2cc(mFn,5);
ccFni    = m2cc(mFni,5);
mmFni    = m2mm(mFni,5);
mmFnit   = m2mm(mFni',5);
mmFnitrc = mmFnit([1 2 3 5 4],[1 2 3 5 4]);
mmFnitr  = mmFnit([1 2 3 5 4],:);
mmFnitc  = mmFnit(:,[1 2 3 5 4]);
mmFBit   = m2mm(mFBi',5);
mmFBitrc = mmFBit([1 2 3 5 4],[1 2 3 5 4]);
mmFt     = m2mm(mF',5);
mmFtr    = mmFt([1 2 3 5 4],:);

mmA  = ccFn * ccGtB';
mmAr = mmA([1 2 3 5 4],:);
mmT  = mmAr + mmA;
mmA1 = 1/2 * ( mmFnitrc + mmFnitc );
mmA2 = mmFni * mmFnitr;
mmP  = mmA1(:,[1 2 3 5 4]) * mmA2;

%======================================================================
% elastic material matrix

mmH  = 2*(Gm - Gl*log(J)) .* mmIs + Gl*ccI*ccI';
mmHc = mmH(:,[1 2 3 5 4]);

%----------------------------------------------------------------------
% The flow criterion is evaluated.

mGtd  = mGt - 1/3*trace(mGt)*mI;
ccGtd = m2cc(mGtd,5);
Gteq  = sqrt( 3/2 * ccGtd' * ccGtd([1 2 3 5 4]) );
F     = Gteq - Gz;

%----------------------------------------------------------------------
if (F>=0)
%----------------------------------------------------------------------
% viscoplastic deformation

cca  = 3/2 * 1/Gteq * ccGtd;
mmb  = 1/Gteq * ( - cca * cca' + 3/2 * mmIs - 1/2 * ccI * ccI');

GfF  = N * (F/Gty)^(N-1) * 1/Gty;
GfF  = 1/2*(GfF + abs(GfF));
FGk  = - h1 - 2*h2*Gk - 3*h3*Gk^2 - 4*h4*Gk^3 - 7*h5*Gk^6;
c1   = (GDt*Gg*GfF)/(1 - GDt*Gg*GfF*FGk);

mmV  = mmI + GDGl .* mmHc * mmb + c1 * mmHc * cca * cca';
mmVi = inv(mmV);

%----------------------------------------------------------------------
else
%----------------------------------------------------------------------
% elastic deformation

GDGl = 0;
mmVi = mmI;
cca  = zeros(5,1);

%----------------------------------------------------------------------
end;
%----------------------------------------------------------------------

mmE  = mmT - 2*Gl * mmIs * (ccen - GDGl * cca) * ccFni' + mmHc * mmP;
mmEr = mmE([1 2 3 5 4],:);
mmC  = 1/J * ( mmVi * mmEr - J * ccGs * ccFni' );

% material matrix for Updated Lagrange formulation

mmM  = mmC * mmFBit * mmFt;  

% stress matrix

mmGst  = m2mm(mGs',5);
mmGstr = mmGst([1 2 3 5 4],:);
mmGS   = ccGs*ccI' - mmGstr;;

%**********************************************************************
