%******************************************************************************

function [mmM,mmGS,ccGs,mF] = plaxelas3(ipda,vrs,mF);

%======================================================================
%  Get integration point values from ipda0.

c0  = ipda(1);
c1  = ipda(2);

%======================================================================
%  Deformation quantities

mI  = eye(3);
mB  = mF * mF';
mA  = (mB - mI)/2;
trA = mA(1,1) + mA(2,2) + mA(3,3);

%======================================================================
%  Plane strain, plane stress, axisymmetric

%----------------------------------------------------------------------
% Note :
% In case of plane stress the deformation tensor is adapted.
% Also the material parameter 'c0' is changed.
% The stress calculation can be done with the original 'c0'(=c00).
% In the consistent stiffness the 'new' 'c0' must be used.
%----------------------------------------------------------------------

if vrs==1,                                   %  plane strain
   mF(3,3) = 1; 
   Gs33 = c0*trA;
   c00 = c0;
elseif vrs==2,                               %  plane stress
   Gs33 = 0;
   cc  = -c0/(c0+c1); 
   mF(3,3) = sqrt(cc*(mB(1,1)+mB(2,2))-2*cc+1); 
   c00 = c0;
   c0 = c0*c1/(c0+c1);
elseif vrs==3,                               %  axisymmetric
   Gs33 = c0*trA + c1*mA(3,3);
   c00 = c0;
end;

mB = mF * mF';

%======================================================================
%  Initialize some variables

ccI   = m2cc(mI,5);
ccIt  = ccI([1 2 3 5 4]);
ccF   = m2cc(mF,5);
mmF   = m2mm(mF,5);
mmFt  = m2mm(mF',5);
mmFtr = mmFt([1 2 3 5 4],:);
mmFrc = mmF([1 2 3 5 4],[1 2 3 5 4]);
mmFc  = mmF(:,[1 2 3 5 4]);
ccB   = m2cc(mB,5);

%======================================================================
%  Stress

ccGs = 1/2 * c00 * ccB' * ccIt * ccI + ...
       1/2 * c1 * ccB - ...
       1/2 * (3*c00 + c1) * ccI;
ccGs(3) = Gs33;

%mGs = c0*trA*mI + c1*mA;
%mGs(3,3) = Gs33;
%ccGs = m2cc(mGs,5);

%======================================================================
%  Consistent material stiffness matrix

mmM = c0 * ccI * ccF' * mmFtr + ...
      1/2 * c1 * ( mmFrc*mmFtr + mmFc*mmFtr );

%======================================================================
%  Stress stiffness matrix

mGs = zeros(3);
mGs(1,1) = ccGs(1); mGs(2,2) = ccGs(2); mGs(3,3) = ccGs(3); 
mGs(1,2) = ccGs(4); mGs(2,1) = ccGs(5);

mmGst  = m2mm(mGs',5);
mmGstr = mmGst([1 2 3 5 4],:);
mmGS   = ccGs*ccI' - mmGstr;

%******************************************************************************
