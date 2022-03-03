%******************************************************************************
%  plaxq8.m
%
%  Generates integration point data for a 8-noded element.
%
%  neip     :  number of integration poins per element (4 or 9)
%  ksi      :  isoparametric coordinates for every integration point
%                 ksi(int.point , 1) = ksi_1
%                 ksi(int.point , 2) = ksi_2
%  psi      :  value of interpolation function in integration points
%                 psi(int.point , node number)
%  psidksi  :  derivative of interpolation functions in integration
%              points
%  ipwf     :  integration point weight factors
%  lokvg    :  location array for degrees of freedom
%
%**********************************************************************

function [ksi,psi,psidksi,ipwf,lokvg] = plaxq8(lok,ne,nndof,neip);

lokvg(1:ne,:) = ... 
  [nndof*(lok(1:ne,3)-1)+1 nndof*(lok(1:ne,3)-1)+2 ...
   nndof*(lok(1:ne,4)-1)+1 nndof*(lok(1:ne,4)-1)+2 ...
   nndof*(lok(1:ne,5)-1)+1 nndof*(lok(1:ne,5)-1)+2 ...
   nndof*(lok(1:ne,6)-1)+1 nndof*(lok(1:ne,6)-1)+2 ...
   nndof*(lok(1:ne,7)-1)+1 nndof*(lok(1:ne,7)-1)+2 ...
   nndof*(lok(1:ne,8)-1)+1 nndof*(lok(1:ne,8)-1)+2 ...
   nndof*(lok(1:ne,9)-1)+1 nndof*(lok(1:ne,9)-1)+2 ...
   nndof*(lok(1:ne,10)-1)+1 nndof*(lok(1:ne,10)-1)+2 ];

if neip==4
  a=sqrt(3)/3;

  ksi(1,1) = -a; ksi(1,2) = -a; ksi(2,1) =  a; ksi(2,2) = -a;
  ksi(3,1) = -a; ksi(3,2) =  a; ksi(4,1) =  a; ksi(4,2) =  a;
elseif neip==9
  a=0.77459;

  ksi(1,1) = -a; ksi(1,2) = -a; ksi(2,1) =  0; ksi(2,2) = -a;
  ksi(3,1) =  a; ksi(3,2) = -a; ksi(4,1) = -a; ksi(4,2) =  0;
  ksi(5,1) =  0; ksi(5,2) =  0; ksi(6,1) =  a; ksi(6,2) =  0;
  ksi(7,1) = -a; ksi(7,2) =  a; ksi(8,1) =  0; ksi(8,2) =  a;
  ksi(9,1) =  a; ksi(9,2) =  a;
end;

for ip=1:neip
  k1 = ksi(ip,1); k2 = ksi(ip,2);

  psi(ip,1) = (k1-1)*(k2-1)*(-k1-k2-1)/4;
  psi(ip,2) = (k1+1)*(k2-1)*(-k1+k2+1)/4; 
  psi(ip,3) = (k1+1)*(k2+1)*( k1+k2-1)/4;
  psi(ip,4) = (k1-1)*(k2+1)*( k1-k2+1)/4;
  psi(ip,5) = (k1*k1-1)*(k2-1)/2;      psi(ip,6) = (-k1-1)*(k2*k2-1)/2;
  psi(ip,7) = (k1*k1-1)*(-k2-1)/2;     psi(ip,8) = (k1-1)*(k2*k2-1)/2;

  k = 2*ip-1; l = k+1;

  psidksi(1,k) = ((k2-1)*(-k1-k2-1) - (k1-1)*(k2-1))/4;
  psidksi(1,l) = ((k1-1)*(-k1-k2-1) - (k1-1)*(k2-1))/4;
  psidksi(2,k) = ((k2-1)*(-k1+k2+1) - (k1+1)*(k2-1))/4;
  psidksi(2,l) = ((k1+1)*(-k1+k2+1) + (k1+1)*(k2-1))/4;
  psidksi(3,k) = ((k2+1)*(k1+k2-1) + (k1+1)*(k2+1))/4;
  psidksi(3,l) = ((k1+1)*(k1+k2-1) + (k1+1)*(k2+1))/4;
  psidksi(4,k) = ((k2+1)*(k1-k2+1) + (k1-1)*(k2+1))/4;
  psidksi(4,l) = ((k1-1)*(k1-k2+1) - (k1-1)*(k2+1))/4;
  psidksi(5,k) = k1*(k2-1);            psidksi(5,l) = (k1*k1-1)/2;
  psidksi(6,k) = -(k2*k2-1)/2;         psidksi(6,l) = (-k1-1)*k2;
  psidksi(7,k) = k1*(-k2-1);           psidksi(7,l) = -(k1*k1-1)/2;
  psidksi(8,k) = (k2*k2-1)/2;          psidksi(8,l) = (k1-1)*k2;
end;

if neip==4
  ipwf = ones(4,1);
elseif neip==9
  p = 0.55556; q = 0.88889;
  ipwf(1) = p*p;        ipwf(2) = p*q;        ipwf(3) = p*p;
  ipwf(4) = p*q;        ipwf(5) = q*q;        ipwf(6) = p*q;
  ipwf(7) = p*p;        ipwf(8) = p*q;        ipwf(9) = p*p;
end;

%******************************************************************************
