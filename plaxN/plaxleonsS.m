%**********************************************************************

function [ipdaC] = plaxleonsS(mF,mFB,GDt,ipda0,ipdaB,ipdaC);

%======================================================================
% Get integration point data from ipda0 and ipdaB

R    = 8.314;
T    = 300;

G    = ipda0(6);
Gk   = ipda0(7);
H    = ipda0(15);

A0   = ipda0(8);
GDH  = ipda0(9);
Gt0  = ipda0(10);
Gm   = ipda0(11);
Dinf = ipda0(13);
h    = ipda0(14);

Gl   = ipdaB(75);
D    = ipdaB(76); 
DB   = ipdaB(76); DJ = DB;
Gh   = ipdaB(77); 
GsVM = ipdaB(78);
GG   = G/Gh;
DD   = h*(1-D/Dinf)*GsVM/(sqrt(6)*Gh);

ctBeB = ipdaB(60:64); mtBeB = zeros(3);
mtBeB(1,1) = ctBeB(1); mtBeB(2,2) = ctBeB(2); mtBeB(3,3) = ctBeB(3); 
mtBeB(1,2) = ctBeB(4); mtBeB(2,1) = ctBeB(5);
mCpnB = inv(mtBeB); mCpnJ = mCpnB;

%======================================================================
% Determine deformation varables

mI   = eye(3);

J    = det(mF);           JB   = det(mFB);
J13  = J^(1/3);           JB13 = JB^(1/3);
mtF  = mF ./ J13;         mtFB = mFB ./ JB13;

mtB  = mtF * mtF';
mtBd = mtB - 1/3*trace(mtB) * mI;

% hydrostatic and hardening stress

msh  = Gk*(J-1) * mI;
mw   = H*mtBd;

% incremental deformation variables

mFn  = mF * inv(mFB);
Jn   = det(mFn);
mtFn = mtF * inv(mtFB);
mtCn = mtFn' * mtFn;

% subincremental variables

[AA,BB] = eig(mtCn);
Gl1 = sqrt(BB(1,1)); Gl2 = sqrt(BB(2,2)); Gl3 = sqrt(BB(3,3));
%fprintf('Gl1 = %9.4g ; Gl2 = %9.4g ; Gl3 = %9.4g \n',Gl1,Gl2,Gl3);
n1 = AA(:,1); n2 = AA(:,2); n3 = AA(:,3);

aa1 = Gl1*Gl1+Gl2*Gl2+Gl3*Gl3;
aa2 = sqrt(1/3*(aa1))-1;
aa3 = aa2/0.00005;
aa4 = ceil(aa3);
nsic = min(max(aa4,30),100);

%fprintf('nsic = %3d ; \r',nsic);
dt = GDt/nsic;
dGl1 = Gl1^(1/nsic); dGl2 = Gl2^(1/nsic); dGl3 = Gl3^(1/nsic); 
dJn = Jn^(1/nsic);

%======================================================================
% Subincremental stress integration
for sic = 1:nsic,
%----------------------------------------------------------------------
%Gl = 1;

Gl1 = dGl1^sic; Gl2 = dGl2^sic; Gl3 = dGl3^sic; 
J   = dJn^sic * JB;
p   = -Gk*(J-1);

mtCn = Gl1*Gl1*n1*n1' + Gl2*Gl2*n2*n2' + Gl3*Gl3*n3*n3';
mtUn = Gl1*n1*n1' + Gl2*n2*n2' + Gl3*n3*n3';

Gdt1 = dt; Gdt2 = dt;

fbibnewtons;                                            % fbibnewtons.m
%newtonSS; D = Dinf + exp(h*GsVM*Gdt1/(Dinf*sqrt(6)*Gh))*(DJ-Dinf);

DJ = D; mCpnJ = mCpn;

%----------------------------------------------------------------------
end;

%mbtBen  = mbtBen ./ (det(mbtBen)^(1/3));
%mbtBend = mbtBen - 1/3*trace(mbtBen)*mI;
%mbsd    = G*mbtBend;

mRn  = mtFn*inv(mtUn);
mtBe = mRn*mbtBen*mRn';
msd  = mRn*mbsd*mRn';
GsVM = sqrt(3/2*trace(msd*msd));       


mGs  = msd + msh + mw;
if Gl>0.999, Gl=1; else, Gl = 1/(1+GDt*GG); end;

%fprintf('%13.10f ',Gl);
%if Gh>A0*exp(GDH/(R*T))*Gt0/1e2, Gh=A0*exp(GDH/(R*T))*Gt0; end;
%fprintf('%13.10f %20.10f \r',Gl,D);

%======================================================================
% Put integration point data in ipdaC

ipdaC(60:64) = [mtBe(1,1) mtBe(2,2) mtBe(3,3) mtBe(1,2) mtBe(2,1)];
ipdaC(65:69) = [mGs(1,1) mGs(2,2) mGs(3,3) mGs(1,2) mGs(2,1)];
ipdaC(70:74) = [mCpn(1,1) mCpn(2,2) mCpn(3,3) mCpn(1,2) mCpn(2,1)];
ipdaC(75)    = Gl;
ipdaC(76)    = D;
ipdaC(77)    = Gh;
ipdaC(78)    = GsVM;

%**********************************************************************
