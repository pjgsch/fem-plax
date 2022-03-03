%**********************************************************************

function [mmM,mmGS,ccGs,mF] = plaxelas41(ipda,vrs,mF);

%======================================================================
%  Get integration point values.

c0  = ipda(1); 
c1  = ipda(2);

%======================================================================
%  Deformation quantities

mI  = eye(3);
mC  = mF' * mF;
mE  = (mC - mI)/2;
trE = mE(1,1) + mE(2,2) + mE(3,3);

%======================================================================
%  Plane strain, plane stress, axisymmetric

if vrs==1                                    %  plane strain
   mF(3,3) = 1;
   me(3,3) = 0;
   P33 = c0*trE;
   c00 = c0;
elseif vrs==2                                %  plane stress
   P33 = 0;
   cc = -c0/(c0+c1);
   mF(3,3) = sqrt(cc*(mC(1,1)+mC(2,2))-2*cc+1); 
   mE(3,3) = 0.5*(mF(3,3)*mF(3,3)-1);
   c00 = c0;
   c0 = c0*c1/(c0+c1);
elseif vrs==3,                               %  axisymmetric
   P33 = c0*trE + c1*mE(3,3);
   c00 = c0;
end;

mC = mF' * mF;
J  = det(mF);
Ji = 1/J;

%======================================================================
%  Initialize some variables

ccI   = m2cc(mI,5);
ccIt  = ccI([1 2 3 5 4]);
ccF   = m2cc(mF,5);
ccFtc = ccF;
mmF   = m2mm(mF,5);
mmFt  = m2mm(mF',5);
mmFtr = mmFt([1 2 3 5 4],:);
mmFrc = mmF([1 2 3 5 4],[1 2 3 5 4]);
%mmFc  = mmF(:,[1 2 3 5 4]);
mmFr  = mmF([1 2 3 5 4],:);
ccC   = m2cc(mC,5);

%======================================================================
%  2nd Piola-Kirchhoff and Cauchy stress 

%ccP = 1/2 * c00 * ccC' * ccIt * ccI + ...
%      1/2 * c1 * ccC - ...
%      1/2 * (3*c00 + c1) * ccI;
%ccP(3) = P33;

mP = c0*trE*mI + c1*mE;
mP(3,3) = P33;
ccP = m2cc(mP,5);

ccPt = ccP([1 2 3 5 4]);
ccGs = Ji * mmF * mmFr * ccPt;

%======================================================================
%  Stress stiffness matrix for Updated Lagrange

mGs = zeros(3);
mGs(1,1) = ccGs(1); mGs(2,2) = ccGs(2); mGs(3,3) = ccGs(3); 
mGs(1,2) = ccGs(4); mGs(2,1) = ccGs(5);

mmGst  = m2mm(mGs',5);
mmGstr = mmGst([1 2 3 5 4],:);
mmGS   = ccGs*ccI' - mmGstr;

%======================================================================
%  Consistent material stiffness matrix for Updated Lagrange

mmM1 = c0 * ccI * ccFtc' * mmFtr + ...
       1/2 * c1 * ( mmFtr*mmFtr + mmFt*mmFtr );
mmGs = m2mm(mGs,5);
mmM  = - ccGs * ccI' + mmGstr + mmGs + Ji * mmF * mmFrc * mmM1;

%**********************************************************************
