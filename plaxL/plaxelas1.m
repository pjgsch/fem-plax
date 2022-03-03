%******************************************************************************
%  Cauchy stresses and linear strains in the 
%  integration points of the elements are calculated.
%  The material behaviour is linear elastic and isotropic.
%----------------------------------------------------------------------
%  function [s,strs,strn,vm] = plaxelas1(eida,vrs,du);
%----------------------------------------------------------------------
%  eida = element integration point data
%  vrs  = 1 : plane strain
%  	= 2 : plane stress
%  	= 3 : axi-symmetric
%  du	= element nodal point displacement gradient
%  s	= matrix with material stiffness matrix (Hooke's law)
%  strs = column with stresses
%  strn = column with strains
%  vm   = Von Mises stress
%**********************************************************************

function [s,strs,strn,vm] = plaxelas1(eida,vrs,du);

em = eida(1);
nu = eida(2);

if vrs==1                                     
  a=em/((1-2*nu)*(1+nu));
  s = a * [
  1-nu   nu     nu     0          0 
  nu     1-nu   nu     0          0
  nu     nu     1-nu   0          0
  0      0      0     (1-2*nu)/2 (1-2*nu)/2
  0      0      0     (1-2*nu)/2 (1-2*nu)/2  ];
elseif vrs==2                                 
  a=em/(1-nu*nu);
  s = a * [
  1   nu  nu  0        0
  nu  1   nu  0        0
  nu  nu  1   0        0
  0   0   0   (1-nu)/2 (1-nu)/2 
  0   0   0   (1-nu)/2 (1-nu)/2  ];
elseif vrs==3                                
  a=em/((1-2*nu)*(1+nu));
  s = a * [
  1-nu    nu     nu     0             0
  nu      1-nu   nu     0             0
  nu      nu     1-nu   0             0
  0       0      0      (1-2*nu)/2    (1-2*nu)/2
  0       0      0      (1-2*nu)/2    (1-2*nu)/2  ];
end;

sp = s * du;

if     vrs==1                          % plane strain
       strs = [sp(1:2) ; sp(3) ; sp(4) ; sp(4)];
       strn = [du(1:2) ; 0 ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5)) ];
elseif vrs==2                          % plane stress
       strs = [sp(1:2) ; 0 ; sp(4) ; sp(4)];
       ez   = -nu/em * (sp(1) + sp(2));
       strn = [du(1:2) ; ez ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5))];
elseif vrs==3                          %  axi-symmetric
       strs = sp;
       strn = [du(1:3) ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5))];
end;

t1 = strs(1).^2 + strs(2).^2 + strs(4).^2;
t2 = strs(1).*strs(2) + strs(1).*strs(4) + strs(2).*strs(4);
t3 = 3*(strs(3).^2);
vm = sqrt(t1-t2+t3);

%******************************************************************************
