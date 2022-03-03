%**********************************************************************
edit plaxL.m

%======================================================================
% Quarter plate with central hole. Axial loading;
% The geometry and loading is specified in 'IPhole.m'.
%----------------------------------------------------------------------
clear all; IPhole; plaxL;
figure; plotmesh([100 1 0 0 0 0 0 0 0],lok,crd0,crd,eidaC,21);
%======================================================================

%======================================================================
% Thick-walled cylinder with internal pressure; planar;
% The coordinates and connectivity are given in 'cyl.m'.
% Type 'edit plaxLcylpl.m' to specify material and loading.
%----------------------------------------------------------------------
clear all; cyl; plaxLcylpl; plaxL;
figure; plotmesh([500 1 0 0 0 0 0 0 0],lok,crd0,crd,[],21);
Sr = crd0(pp(1:9,1),2); Sur = Mp(pp(1:9,1),2); 
figure; plot(Sr,Sur);grid on;xlabel('r [m]');ylabel('u_r [m]');
%======================================================================

%======================================================================
% Thick-walled cylinder with internal pressure; axisymmetric
% Type 'edit plaxLcylax.m' to specify matrial and loading.
%----------------------------------------------------------------------
clear all; plaxLcylax; plaxL; 
for e=1:ne
  Srm(e) = (1/2)*( crd0(lok(e,3),1) + crd0(lok(e,4),1) );
  SGsrr(e) = (1/2)* ( eidaC(4*(e-1)+1,21) + eidaC(4*(e-1)+2,21) );
  SGstt(e) = (1/2)* ( eidaC(4*(e-1)+1,23) + eidaC(4*(e-1)+2,23) );
  SGszz(e) = (1/2)* ( eidaC(4*(e-1)+1,22) + eidaC(4*(e-1)+2,22) );
end;
for e=1:ne
  Sr(e) = crd0(lok(e,3),1); Sur(e) = Mp(lok(e,3),1);
end;

figure; plotmesh([1000 1 0 0 0 0 0 0 0],lok,crd0,crd,[],21);

Sr(e+1) = crd0(lok(e,4),1); Sur(e+1) = Mp(lok(e,4),1);

figure; plot(Sr,Sur);grid on;xlabel('r [m]');ylabel('u_r [m]');
figure; plot(Srm,SGsrr,'-k',Srm,SGstt,'--k',Srm,SGszz,'-.k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');ax=axis;axis([0 b ax(3) ax(4)]);
        legend('\sigma_{rr}','\sigma_{tt}','\sigma_{zz}');                              
%======================================================================

%======================================================================
% Thick-walled cylinder with internal pressure; axisymmetric
% Type 'edit IPcylpi.m' to specify geometry, material and loading.
%----------------------------------------------------------------------
clear all; IPcylpi; plaxL;

figure; plotmesh([1000 1 0 0 0 0 0 0 0],lok,crd0,crd,eidaC,21);

for i=1:ne, xxx(2*(i-1)+1)=4*(i-1)+1; xxx(2*i)=4*(i-1)+2; end;
 SGsrr = eidaC(xxx,21); SGszz = eidaC(xxx,22);
 SGstt = eidaC(xxx,23); Sr    = eidaC(xxx,90);

figure; plot(Sr,SGsrr,Sr,SGszz,Sr,SGstt);grid on;xlabel('r');ylabel('\sigma'); 
        legend('\sigma_{rr}','\sigma_{zz}','\sigma_{tt}');
%======================================================================

%======================================================================
% Tensile; linear elastic; plane stress                 edit plaxLtns.m
%----------------------------------------------------------------------
clear all; ety=11; plaxLtns; plaxL;
Mp, Mfe, Mfi, eidaC(1,22), eidaC(1,18)
figure; plotmesh([10 1 1 0 0 0 0 0 0],lok,crd0,crd,[],2);
%----------------------------------------------------------------------
% Tensile; linear elastic; axisymmetric                 edit plaxLtns.m
%----------------------------------------------------------------------
clear all; ety=10; plaxLtns; plaxL;
Mp, Mfe, Mfi, eidaC(1,22), eidaC(1,18)
figure; plotmesh([10 1 1 0 0 0 0 0 0],lok,crd0,crd,[],2);
%======================================================================
% Simple shear; linear elastic;                         edit plaxLshr.m
%----------------------------------------------------------------------
clear all; ety=11; plaxLshr; plaxL;
Mp, Mfe, Mfi, eidaC(1,24), eidaC(1,20)
figure; plotmesh([10 1 1 0 0 0 0 0 0],lok,crd0,crd,[],2);
%======================================================================

%**********************************************************************
