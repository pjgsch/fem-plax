%**********************************************************************
%loka  = lok;
nel   = size(loka,1);  
nenod = size(loka,2) - 2;  
nndof = 2;
neip  = 4;
[ksi,psi,psidksi,ipwf,lokvg] = plaxq4(loka,nel,nndof,neip);

for e=1:nel
  ec0 = crd0(loka(e,3:nenod+2),:);
  ec  = crd(loka(e,3:nenod+2),:); 
  eidaC(neip*(e-1)+1:neip*(e-1)+neip,26:27) = psi*ec0;
  eidaC(neip*(e-1)+1:neip*(e-1)+neip,28:29) = psi*ec;

  ipcrd = eidaC(neip*(e-1)+1:neip*e,28:29);
  n111 = ipcrd(2,1) - ipcrd(1,1);  n112 = ipcrd(2,2) - ipcrd(1,2);
  l11  = sqrt(n111*n111+n112*n112);
  n111 = n111/l11;  n112 = n112/l11;
  n121 = ipcrd(3,1) - ipcrd(1,1);  n122 = ipcrd(3,2) - ipcrd(1,2);
  l12  = sqrt(n121*n121+n122*n122);
  n121 = n121/l12;  n122 = n122/l12;
  n221 = ipcrd(4,1) - ipcrd(2,1);  n222 = ipcrd(4,2) - ipcrd(2,2);
  l22  = sqrt(n221*n221+n222*n222);
  n221 = n221/l22;  n222 = n222/l22;
  n311 = ipcrd(4,1) - ipcrd(3,1);  n312 = ipcrd(4,2) - ipcrd(3,2);
  l31  = sqrt(n311*n311+n312*n312);
  n311 = n311/l31;  n312 = n312/l31;
  eidaC(neip*(e-1)+1:1:neip*e,30:33) = ...
  [ n111 n112 n121 n122 ; n111 n112 n221 n222
    n311 n312 n121 n122 ; n311 n312 n221 n222 ];
end;



%**********************************************************************
