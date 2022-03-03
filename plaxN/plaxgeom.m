%**********************************************************************

function  [thk0,thkB,thk,dt0,dtB,dt,dfie0,dfieB,dfie,mF,mFB,r0,r] = ...
          plaxgeom(ec0,ecB,ec,psi,psidksi,vrs,ipda0,ip,nenod);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates geometrical quantities for in integration point.
% It also calculates the deformation matrix.
% Values are calculated in the initial state (0), the begin-increment 
% state (B) and the current state (-).
% A difference is made between plane strain/stress (vrs=1|2) and 
% axisymmetry (vrs=3).
%
% output
%  thk     = thickness
%  dt      = determinant of the Jacobian matrix
%  dfie    = beta-matrix : derivatives of shape functions
%  mF      = deformation matrix
%  r       = radius (in case of axisymmetry)
% input
%  ec      = nodal point coordinates
%  psi     = shape functions 
%  psidksi = derivatives of shape functions w.r.t. local coordinates
%  vrs     = plane strain (1), plane stress (2), axisymm. (3)
%  ipda0   = initial integration point data
%  ip      = integration point number
%  nenod   = number of element nodes
% local  
%  ax      = either 0 or 1
%  dpsi    = derivatives of shape functions w.r.t. local coordinates
%  dpsixy  = derivatives of shape functions w.r.t. global coordinates
%  jc      = Jacobian matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if vrs==3
     r0      = psi(ip,:)*ec0(:,1);
     rB      = psi(ip,:)*ecB(:,1);
     r       = psi(ip,:)*ec(:,1);
     thk0    = r0*2*pi;
     thkB    = rB*2*pi;
     thk     = r*2*pi;
     mF33    = r/r0;
     mFB33   = rB/r0;
     ax      = 1;
  else
     r0      = 1;             rB      = 1;        r       = 1;
     thk0    = ipda0(3);      thkB    = thk0;     thk     = thk0;
     mF33    = 1;             mFB33   = 1;
     ax      = 0;
  end;

  dpsi(:,1) = psidksi(:,2*ip-1); 
  dpsi(:,2) = psidksi(:,2*ip);

  jc0  = dpsi' * ec0;	  jci0 = inv(jc0);
  jcB  = dpsi' * ecB;     jciB = inv(jcB);
  jc   = dpsi' * ec;      jci  = inv(jc);
  dt   = det(jc);         dt0  = det(jc0);      dtB  = det(jcB);   
  
  dfie   = zeros(5,2*nenod);
  dpsixy = dpsi * jci' ;
  dfie(1,2*(1:nenod)-1) = dpsixy(1:nenod,1)'; 
  dfie(2,2*(1:nenod))	= dpsixy(1:nenod,2)'; 
  dfie(3,2*(1:nenod)-1) = ax.*psi(ip,1:nenod)/r;
  dfie(4,2*(1:nenod))	= dpsixy(1:nenod,1)';
  dfie(5,2*(1:nenod)-1) = dpsixy(1:nenod,2)'; 

  dfie0   = zeros(5,2*nenod);
  dpsixy0 = dpsi * jci0' ;
  dfie0(1,2*(1:nenod)-1) = dpsixy0(1:nenod,1)'; 
  dfie0(2,2*(1:nenod))   = dpsixy0(1:nenod,2)'; 
  dfie0(3,2*(1:nenod)-1) = ax.*psi(ip,1:nenod)/r0;
  dfie0(4,2*(1:nenod))   = dpsixy0(1:nenod,1)';
  dfie0(5,2*(1:nenod)-1) = dpsixy0(1:nenod,2)'; 

  dfieB   = zeros(5,2*nenod);
  dpsixyB = dpsi * jciB' ;
  dfieB(1,2*(1:nenod)-1) = dpsixyB(1:nenod,1)'; 
  dfieB(2,2*(1:nenod))   = dpsixyB(1:nenod,2)'; 
  dfieB(3,2*(1:nenod)-1) = ax.*psi(ip,1:nenod)/rB;
  dfieB(4,2*(1:nenod))   = dpsixyB(1:nenod,1)';
  dfieB(5,2*(1:nenod)-1) = dpsixyB(1:nenod,2)'; 

  mF  = eye(3); mF(1:2,1:2)  = jci0*jc;  mF(3,3)  = mF33;  mF  = mF';
  mFB = eye(3); mFB(1:2,1:2) = jci0*jcB; mFB(3,3) = mFB33; mFB = mFB';

%**********************************************************************
