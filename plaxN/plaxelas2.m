%******************************************************************************
%  Material parameters for orthotropic linear elastic behaviour are
%  extracted from the database array 'eida'.
%  Unknown Poisson ratio's are calculated from the symmetry conditions.
%----------------------------------------------------------------------
%  function [sg,sl,ssl,ssg,snl,sng] = plaxelas2(eida,vrs,du);
%----------------------------------------------------------------------
%  eida = element integration point data
%  vrs  = 1 : plane strain
%  	= 2 : plane stress
%  	= 3 : axi-symmetric
%  du	= element nodal point displacement gradient
%**********************************************************************

function [sg,sl,ssl,ssg,snl,sng] = plaxelas2(eida,vrs,du);

e11 = eida(1);  
n12 = eida(2);  
e22 = eida(3);  
g12 = eida(4);  
e33 = eida(5);  
n31 = eida(6);  
n32 = eida(7);  
ang = eida(8);  

n21 = n12*e22/e11; 
%n12 = n21*e11/e22; 
n13 = n31*e11/e33; 
%n13 = n31*e11/e33; 
n23 = n32*e22/e33;
%n32 = n23*e33/e22;

% Some constants are calculated.

a   = 1-n21*n12;
d   = e11*e22*e33/(1-n12*n21-n23*n32-n31*n13-n12*n23*n31-n21*n32*n13);

% The material stiffness matrix in the local, matarial coordinate
% system is 'sl'.
% It is first initialised and then made for plane strain (vrs=1),
% plane stress (vrs=2) and axisymmetry (vrs=3).

sl = zeros(5);

if vrs==1                              % plane strain
  sl(1,1) = d*(1-n32*n23)/(e22*e33);
  sl(1,2) = d*(n31*n23+n21)/(e22*e33);
  sl(1,3) = d*(n21*n32+n31)/(e22*e33);
%  sl(1,3) = 0;
  sl(2,1) = d*(n13*n32+n12)/(e11*e33);
  sl(2,2) = d*(1-n31*n13)/(e11*e33);
  sl(2,3) = d*(n12*n31+n32)/(e11*e33);
%  sl(2,3) = 0;
  sl(3,1) = d*(n12*n23+n13)/(e11*e22);
  sl(3,2) = d*(n21*n13+n23)/(e11*e22);
  sl(3,3) = d*(1-n12*n21)/(e11*e22);
  sl(4,4) = g12;
  sl(4,5) = g12;
  sl(5,4) = g12;
  sl(5,5) = g12;
elseif vrs==2                          % plane stress
  sl(1,1) = e11/a;
  sl(1,2) = n21*e11/a;
  sl(1,3) = 0; 
  sl(2,1) = n12*e22/a; 
  sl(2,2) = e22/a;   
  sl(2,3) = 0;
  sl(3,1) = 0;
  sl(3,2) = 0;
  sl(3,3) = 0;
  sl(4,4) = g12;
  sl(4,5) = g12;
  sl(5,4) = g12;
  sl(5,5) = g12;
elseif vrs==3                          % axi-symmetric
  sl(1,1) = d*(1-n32*n23)/(e22*e33);
  sl(1,2) = d*(n31*n23+n21)/(e22*e33);
  sl(1,3) = d*(n21*n32+n31)/(e22*e33);
  sl(2,1) = d*(n13*n32+v12)/(e11*e33);
  sl(2,2) = d*(1-n31*n13)/(e11*e33);
  sl(2,3) = d*(n12*n31+n32)/(e11*e33);
  sl(3,1) = d*(n12*n23+n13)/(e11*e22);
  sl(3,2) = d*(n21*n13+n23)/(e11*e22);
  sl(3,3) = d*(1-n12*n21)/(e11*e22);
  sl(4,4) = d*g12;
  sl(4,5) = d*g12;
  sl(4,5) = d*g12;
  sl(5,5) = d*g12;
end;

% To calculate the material stiffness matrix in the global coordinate
% system ('sg') the local matrix 'sl' must be transformed.
% The transformation matrix 'tm' is made using cosine and sine of the 
% angle 'ang' between the local '1'-axis and the global 'x'-axis.

co = cos(ang*pi/180); si = sin(ang*pi/180);

tm = [ ...
co*co    si*si   0 -co*si      -co*si
si*si    co*co   0  co*si       co*si
0        0       1  0           0
co*si   -co*si   0  co*co      -si*si 
co*si   -co*si   0 -si*si       co*co ];

% Transformation is carried out and therefor the inverse of 'tm'
% is needed.

tmi = inv(tm);

sg = tmi' * sl * tmi;

%======================================================================
% Now stresses are calculated.
%
% Global stress components 'sp' are simply calculated by multiplying 
% the global material stiffness matrix 'sg' with the derivatives of the
% displacements 'du'.

sp = sg * du;

% The stresses and strains are stored in columns for 
% plane strain (vrs=1), plane stress (vrs=2) and axisymmetry (vrs=3).

if     vrs==1                          % plane strain

       ssg = [sp(1:2) ; sp(3) ; sp(4) ; sp(4)];
       sng = [du(1:2) ; 0 ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5))];

elseif vrs==2                          % plane stress

       ssg = [sp(1:2) ; 0 ; sp(4) ; sp(4)];
       ez  = -n12/e11*sp(1) - n23/e22*sp(2);
       sng = [du(1:2) ; ez ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5))];

elseif vrs==3                          %  axi-symmetric

       ssg = sp;
       sng = [du(1:3) ; 0.5*(du(4) + du(5)) ; 0.5*(du(4) + du(5))];

end;

% Stress and strain components in the material coordinate system
% 'ssl' and 'snl' are calculated by transformation using the inverse 
% transformation matrix 'tmi'.

ssl = tmi * ssg;
snl = tmi * sng;

%******************************************************************************
