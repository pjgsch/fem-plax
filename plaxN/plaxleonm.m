%**********************************************************************

function [mmM,mmGS] = plaxleonm(mF,mFB,GDt,ipda0,ipdaB,ipdaC);

%======================================================================
%  Get integration point values from ipda0, ipdaC and ipdaC

mI   = eye(3);
ccI  = m2cc(mI,5);
mmI  = m2mm(mI,5);

G    = ipda0(6);
Gk   = ipda0(7);
H    = ipda0(15);

A0   = ipda0(8);
GDH  = ipda0(9);
Gt0  = ipda0(10);
Gm   = ipda0(11);
Dinf = ipda0(13);
h    = ipda0(14);

Gl   = ipdaC(75);
D    = ipdaC(76);
Gh   = ipdaC(77);
GsVM = ipdaC(78);

ctBeB = ipdaB(60:64); mtBeB = zeros(3);
mtBeB(1,1) = ctBeB(1); mtBeB(2,2) = ctBeB(2); mtBeB(3,3) = ctBeB(3); 
mtBeB(1,2) = ctBeB(4); mtBeB(2,1) = ctBeB(5);
%mtbBeiB = inv(mtBeB);

ctBe = ipdaC(60:64); mtBe = zeros(3);
mtBe(1,1) = ctBe(1); mtBe(2,2) = ctBe(2); mtBe(3,3) = ctBe(3); 
mtBe(1,2) = ctBe(4); mtBe(2,1) = ctBe(5);
mtBed = mtBe -1/3*trace(mtBe)*mI;

ccGs = ipdaC(65:69)'; mGs = zeros(3);
mGs(1,1) = ccGs(1); mGs(2,2) = ccGs(2); mGs(3,3) = ccGs(3); 
mGs(1,2) = ccGs(4); mGs(2,1) = ccGs(5);

%======================================================================
%  Initialize some variables

J   = det(mF);       JB   = det(mFB);
J13 = J^(1/3);       JB13 = JB^(1/3);
mtF = mF ./ J13;     mtFB = mFB ./ JB13;
mtC = mtF'*mtF;      mCpB = mtFB'*inv(mtBeB)*mtFB;

mCp     = (1-Gl)*mtC + Gl*mCpB;
mFi     = inv(mF);
mCpi    = inv(mCp);

ccF     = m2cc(mF,5);
ccFi    = m2cc(mFi,5);
ccFit   = ccFi([1 2 3 5 4]);
mmtF    = m2mm(mtF,5);
mmtFcr  = mmtF([1 2 3 5 4],[1 2 3 5 4]);
mmtFt   = m2mm(mtF',5);
mmtFtr  = mmtFt([1 2 3 5 4],:);
mmtFc   = mmtF(:,[1 2 3 5 4]);
mmFt    = m2mm(mF',5);
mmFBt   = m2mm(mFB',5);
ccCp    = m2cc(mCp,5);
ccCpB   = m2cc(mCpB,5);
cctC    = m2cc(mtC,5);
cctBedt = m2cc(mtBed',5);

%======================================================================
%  Variation of strain variables

%  variation of isochoric elastic left Cauchy strain tensor

mM1    = mtF * mCpi';
mmM1   = m2mm(mM1,5);
mmM1c  = mmM1(:,[1 2 3 5 4]);
mmM1cr = mmM1c([1 2 3 5 4],:);
mM2    = mtF * mCpi;
mmM2   = m2mm(mM2,5);
mmM2c  = mmM2(:,[1 2 3 5 4]);

mmA1   = mmM1cr + mmM2c;
mmA2   = - mmM2c * mmM1c;

%  variation of isochoric deformation tensor

mmF    = -(1/3)/J13 * ccF * ccFit' + 1/J13 * mmI;

%  variation of plastic right Cauchy strain

mmC1   = (1-Gl)*(mmtFtr + mmtFt);
ccC2   = (ccCpB - cctC);

%  variation of elastic parameter

if GsVM>0
ce11 = exp(GsVM/(sqrt(3)*Gt0));
ce12 = exp(-GsVM/(sqrt(3)*Gt0));
%ce1  = Gh*(1/GsVM - 1/(sqrt(3)*Gt0));
ce1  = Gh/GsVM + Gh*(ce11+ce12)/(ce11-ce12);

DD   = h*(1-D/Dinf)*GsVM/(sqrt(6)*Gh);    
%cdd1 = DD/(sqrt(3)*Gt0);
cdd1 = DD*(1/GsVM - ce1/Gh);
cdd2 = -DD*Gm/Gt0;
cdd3 = DD - h*GsVM/(sqrt(6)*Dinf*Gh);
cd1  = GDt*cdd1/(1-GDt*cdd3);
cd2  = GDt*cdd2/(1-GDt*cdd3);
ce2  = Gh*Gm/Gt0;
ce3  = -Gh;
GG   = G/Gh;
ll   = Gl*GDt*GG/(GDt*G + Gh);
h1   = 3*G*G/(2*GsVM) * (ce1 + ce3*cd1);
h2   = -Gk*J* (ce2 + ce3*cd2);
l1   = ll*h1;
l2   = ll*h2;
else
l1   = 0;
l2   = 0;
end;

%l2=0;
%======================================================================
% Deviatoric stress

mmB1   = mmA1 + mmA2*mmC1;
ccB2   = mmA2 * ccC2;

mmH1   = mmI - 1/3*ccI*ccI';
mmH2   = inv(mmI - l1*ccB2*cctBedt');
mmH3   = mmB1*mmF + l2*ccB2*ccI';

mmSd   = G*mmH1*mmH2*mmH3;

%======================================================================
% Hydrostatic stress 

mmSh   = Gk*ccI*J*(ccFit');

%======================================================================
% Hardening stress 

mmH    = H*(mmtFcr + mmtFc - 2/3*ccI*ccI'*mmtFc)*mmF;

%======================================================================
% Consistent material stiffness matrix

mmS    = mmSd + mmSh + mmH;
mmSc   = mmS(:,[1 2 3 5 4]);
mmM    = mmSc*mmFt;

%======================================================================
% Stress stiffness matrix

mmGst  = m2mm(mGs',5);
mmGstr = mmGst([1 2 3 5 4],:);
mmGS   = ccGs*ccI' - mmGstr;

%**********************************************************************
