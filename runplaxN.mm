%**********************************************************************
edit plax.m

%======================================================================
% Tensile test; linear elastic;                          
% Type 'edit plaxtns.m' to specify material, geometry and loading.
%----------------------------------------------------------------------
clear all; ety=11;ma=1;A0=1; nic=200; plaxtns; plax; plaxpostdata;
figure; plotmesh([10 1 1 0 0 0 0 0 0],lok,crd0,crd,[],2);
%======================================================================

%======================================================================
% Shear test; linear elastic;
% Type 'edit plaxshr.m' to specify input.
%----------------------------------------------------------------------
clear all; ety=11;ma=1;A0=1; plaxshr; plax; plaxpostdata;
figure; plotmesh([10 1 1 0 0 0 0 0 0],lok,crd0,crd,[],2);
%======================================================================

%======================================================================
% Tensile test; linear elastic; plane stress      
% Type 'edit plaxGsGeltnsps.m' to specify input.
%----------------------------------------------------------------------
clear all; plaxGsGeltnsps; plax; plaxpostdata; 
figure; plotmesh([20 1 1 0 0 0 0 0 0],lok,crd0,crd,eidaC,2);
%======================================================================

%======================================================================
% Inhomogeneous deformation 
% Type 'edit plaxinhdef.m' to see the input.
% Deformation is prescribed and the deformation tensor is calculated.
% The engineering stress is calculated with material model mat=1, which 
% is a linear relation between engineering stress and linear strain.
%----------------------------------------------------------------------
clear all; ety=11; ma=1; plaxinhdef; plax; plaxpostdata;                               
figure; plotmesh([1 1 1 0 0 0 0 0 0],lok,crd0,crd,eidaC,2);
%----------------------------------------------------------------------
% Analytical solution; deformation matrix is printed;
%----------------------------------------------------------------------
l0=crd0(2,1)-crd0(1,1); h0=crd0(3,2)-crd0(2,2); 
ll=crd(2,1)-crd(1,1); hh=crd(3,2)-crd(2,2);
for i=1:4, 
 x01=eidaB(i,26); x02=eidaB(i,27);
 F11=(ll/l0); F12=0; F21=(hh-h0)/(h0*l0)*x02; F22=1+(hh-h0)/(h0*l0)*x01;
 F33=1; MF=[F11 F12 0; F21 F22 0; 0 0 1]; deF=det(MF);
 SF(i,1)=F11;SF(i,2)=F22;SF(i,3)=F33;SF(i,4)=F12;SF(i,5)=F21;SF(i,6)=deF;
end;
for i=1:4, fprintf('%11.4e %11.4e %11.4e %11.4e %11.4e %11.4e \n',SF(i,1:6));end;
%======================================================================

%======================================================================
% Inhomogeneous deformation.
% Type 'edit plaxinhdef.m' to see the input.
% Deformation is prescribed and the deformation tensor is calculated.
% The Cauchy stress is calculated with material model mat=3,
% which is a linear relation between Cauchy stress and Finger strain.
%----------------------------------------------------------------------
clear all; ety=11; ma=3; plaxinhdef; plax; plaxpostdata;
figure; plotmesh([1 1 1 0 0 0 0 0 0],lok,crd0,crd,eidaC,2);
%----------------------------------------------------------------------
% Analytical solution; Cauchy stress matrix is printed.
%----------------------------------------------------------------------
l0=crd0(2,1)-crd0(1,1); h0=crd0(3,2)-crd0(2,2); 
ll=crd(2,1)-crd(1,1); hh=crd(3,2)-crd(2,2);
C=elda(1,6); Gn=elda(1,7); c0=C*Gn/((1-2*Gn)*(1+Gn)); c1=C/(1+Gn);
for i=1:4, 
 x01=eidaB(i,26); x02=eidaB(i,27);
 F11=(ll/l0); F12=0; F21=(hh-h0)/(h0*l0)*x02; F22=1+(hh-h0)/(h0*l0)*x01;
 F33=1; MF=[F11 F12 0; F21 F22 0; 0 0 1]; deF=det(MF);
 mI=eye(3); mB=MF*MF'; mA=(mB-mI)/2; trA=mA(1,1)+mA(2,2)+mA(3,3); 
 MS=c0*trA*mI+c1*mA;
 SS(i,1)=MS(1,1);SS(i,2)=MS(2,2);SS(i,3)=MS(3,3);
 SS(i,4)=MS(1,2);SS(i,5)=MS(2,1);
end;
for i=1:4, fprintf('%8.2e %8.2e %8.2e %8.2e \n',SS(i,1:4));end;
%======================================================================

%======================================================================
% Inhomogeneous deformation.
% Type 'edit plaxinhdef.m' to see the input.
% Deformation is prescribed and the deformation tensor is calculated.
% The 2nd PK stress is calculated with material model mat=4, which 
% is a linear relation between 2ndPK stress and Green-Lagrange strain.
% The Cauchy stress is calculated from the 2nd PK stress.
%----------------------------------------------------------------------
clear all; ety=11; ma=4; plaxinhdef; plax; plaxpostdata;
figure; plotmesh([1 1 1 0 0 0 0 0 0],lok,crd0,crd,eidaC,2);
%----------------------------------------------------------------------
% Analytical solution; 
% 2nd Piola-Kirchhoff and Cauchy stress matrices are printed.
%----------------------------------------------------------------------
l0=crd0(2,1)-crd0(1,1); h0=crd0(3,2)-crd0(2,2); 
ll=crd(2,1)-crd(1,1); hh=crd(3,2)-crd(2,2);
C=elda(1,6); Gn=elda(1,7); c0=C*Gn/((1-2*Gn)*(1+Gn)); c1=C/(1+Gn);
for i=1:4, 
 x01=eidaB(i,26); x02=eidaB(i,27);
 F11=(ll/l0); F12=0; F21=(hh-h0)/(h0*l0)*x02; F22=1+(hh-h0)/(h0*l0)*x01;
 F33=1; MF=[F11 F12 0; F21 F22 0; 0 0 1]; deF=det(MF);
 mI=eye(3); mC=MF'*MF; mE=(mC-mI)/2; trE=mE(1,1)+mE(2,2)+mE(3,3); 
 MP=c0*trE*mI+c1*mE; MS=(1/deF)*MF*MP*MF';
 SP(i,1)=MP(1,1);SP(i,2)=MP(2,2);SP(i,3)=MP(3,3);
 SP(i,4)=MP(1,2);SP(i,5)=MP(2,1);
 SS(i,1)=MS(1,1);SS(i,2)=MS(2,2);SS(i,3)=MS(3,3);
 SS(i,4)=MS(1,2);SS(i,5)=MS(2,1);
end;
for i=1:4, fprintf('%11.4e %11.4e %11.4e %11.4e \n',SP(i,1:4)); end;
for i=1:4, fprintf('%11.4e %11.4e %11.4e %11.4e \n',SS(i,1:4)); end;
%======================================================================

%======================================================================
% Tensile test; elastic; Select proper model;        edit plaxelastns.m
%----------------------------------------------------------------------
clear all; ma=1;ety=3; A0=0.1;nic=250;tol=0;plaxelatns; plax;
clear all; ma=3;ety=3; A0=0.1;nic=150;tol=0;plaxelatns; plax;
clear all; ma=3;ety=11;A0=0.1;nic=133;tol=0;plaxelatns; plax;
clear all; ma=3;ety=3; A0=0.1;nic=175;tol=1;plaxelatns; plax;
clear all; ma=4;ety=3; A0=0.1;nic=150;tol=0;plaxelatns; plax;
clear all; ma=4;ety=3; A0=0.1;nic=150;tol=1;plaxelatns; plax;
clear all; ma=4;ety=11;A0=0.1;nic=150;tol=0;plaxelatns; plax;
%----------------------------------------------------------------------
figure; plotplot(SGl,Sfy,2);grid on;
        xlabel('\lambda');ylabel('F [N]');
figure; plotplot(SGl,SGs22,2);grid on;
        xlabel('\lambda');ylabel('\sigma_{22} [MPa]');        
figure; plotplot(SGl,SA,2);grid on; 
        xlabel('\lambda');ylabel('A [mm^2]');
figure; mshop = [1 1 0 0 0 0 0 0 0]; ax=[-1 200 -1 300]; plotmovie;
%======================================================================

%======================================================================
% Rigid rotation of one element;                      edit plaxrigrot.m
%----------------------------------------------------------------------
plaxrigrot; plax; 
%----------------------------------------------------------------------
figure; plot(Sfi,Scx); grid on; xlabel('\phi'); ylabel('crd-1 node 3');
figure; plot(Sfi,Sff); grid on; xlabel('\phi'); ylabel('force node 3');
figure; ax=[-2 2 -2 2]; mshop = [1 1 0 0 0 0 0 0 0]; ps=0.01; plotmovie;
%======================================================================

%======================================================================
% Simple shear; elastic;                              edit plaxelashr.m
%----------------------------------------------------------------------
clear all; ma=3;ety=3; A0=0.1;nic=250;tol=0;plaxelashr; plax;     
clear all; ma=3;ety=11;A0=0.1;nic=250;tol=0;plaxelashr; plax;
clear all; ma=4;ety=3; A0=0.1;nic=250;tol=0;plaxelashr; plax; 
clear all; ma=4;ety=3; A0=0.1;nic=250;tol=1;plaxelashr; plax;
clear all; ma=4;ety=11;A0=0.1;nic=250;tol=0;plaxelashr; plax;
%----------------------------------------------------------------------
figure; plot(SGg,Sfiy);grid on;xlabel('\gamma');ylabel('F_y [N]');      
figure; plot(SGg,Sfix);grid on;xlabel('\gamma');ylabel('F_x [N]');      
figure; ax=[-100 400 -1 200]; plotmovie;
%======================================================================

%======================================================================
% Tensile test; linear viscoelastic; multimode        edit plaxvietns.m
%----------------------------------------------------------------------
clear all; ma=8;ety=10;A0=10;nic=99;GDt=0.01;plaxvietns; plax;
%----------------------------------------------------------------------
figure; plot(Sti,SGs22,'b');grid on;
        xlabel('t [s]');ylabel('\sigma [MPa]'); 
%======================================================================

%======================================================================
% Tensile test; Perzyna;                              edit plaxprztns.m
%----------------------------------------------------------------------
clear all; ma=5;ety=3;A0=0.1;nic=100;GDt=0.001;Gdu=1;plaxprztns; plax;
clear all; ma=5;ety=3;A0=0.1;nic=100;GDt=0.01; Gdu=1;plaxprztns; plax;
clear all; ma=5;ety=3;A0=0.1;nic=100;GDt=0.1;  Gdu=1;plaxprztns; plax;
%----------------------------------------------------------------------
SGl = (crd0(3,2)+Su32)/crd0(3,2); 
Sfy = Sfi32+Sfi42; 
de = (SGl(20)-SGl(19))/GDt; 
figure; plot(SGl,S1161);grid on;xlabel('\lambda');ylabel('\sigma [MPa]');
figure; plot(SGl,Sfy);grid on;xlabel('\lambda');ylabel('F [N]'); 
figure; mshop = [1 1 1 0 0 0 0 0 0]; ax=[-1 100 -1 300]; plotmovie;
%======================================================================

%======================================================================
% Simple shear test; Perzyna;                         edit plaxprzshr.m
%----------------------------------------------------------------------
clear all; ma=5;ety=11;A0=1;nic=100;GDt=1;Gdu=1;plaxprzshr; plax; 
%----------------------------------------------------------------------
de = ( STp31(20)/crd0(3,2) - STp31(19)/crd0(3,2) )/GDt;
figure; plotplot(STp31/crd0(3,2),(Sfi31+Sfi41),2);grid on;
        xlabel('\gamma');ylabel('F_x [N]');
figure; mshop = [1 1 1 0 0 0 0 0 0];plotmesh(mshop,lok,crd0,crd,eidaC,21);
figure; mshop = [1 1 0 0 0 0 0 0 0]; ax=[-50 200 -10 100]; plotmovie;
%======================================================================

%======================================================================
% Tensile test; EGP;                                  edit plaxegptns.m
%----------------------------------------------------------------------
clear all; ma=6;ety=11;A0=0.1;nic=400;GDt=0.1;Gdu=1;plaxegptns; plax;
%----------------------------------------------------------------------
SGl = (crd0(3,2)+Su32)/crd0(3,2); 
Sfy = Sfi32+Sfi42; 
de = (SGl(20)-SGl(19))/GDt;
figure; plot(SGl,S1161);grid on;xlabel('\lambda');ylabel('\sigma [MPa]');
figure; plot(SGl,Sfy);grid on;xlabel('\lambda');ylabel('F [N]');
figure; mshop = [1 1 1 0 0 0 0 0 0]; ax=[-1 100 -1 300];  plotmovie;
%======================================================================

%======================================================================
% Simple shear test; EGP;                             edit plaxegpshr.m
%----------------------------------------------------------------------
clear all; ma=6;ety=11;A0=0.1;nic=100;GDt=1;Gdu=1;plaxegpshr; plax;
%----------------------------------------------------------------------
SGl = STp31/crd0(3,2); Sfx = (Sfi31+Sfi41); 
figure; plot(SGl,S1161);grid on;xlabel('\lambda');ylabel('\sigma [MPa]');
figure; plot(SGl,Sfx);grid on;xlabel('\lambda');ylabel('F [N]'); 
figure; plotmesh([1 1 1 0 0 0 0 0 0],lok,crd0,crd,eidaC,21);axis([-50 200 -10 100]);
figure; mshop = [1 1 1 0 0 0 0 0 0]; ax=[-1 100 -1 300];  plotmovie;
%======================================================================

%**********************************************************************
