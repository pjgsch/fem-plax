%**********************************************************************
%path('../femaxi',path);
edit femaxi.m

%======================================================================
% Disk with central hole, isotropic, plane stress, edge displacement
%----------------------------------------------------------------------
close all; clear all; 
sta = 'pss'; ne = 50; elt = 2;
a = 0.25; b = 0.5; Gr = 7500; E = 250e9; Gn = 0.33;
Go = 0; Ga = 0;
pf = [ ]; pu = [ ne+1 0.01 ];
Ri = a; Ro = b;
% Element data base 
%     [ E1  E2   E3   Gn12   Gn21   Gn23   Gn32   Gn31  Gn13   G12 ]
EE1 = [ E   E    E    Gn     0      0      Gn     Gn    0      E/(2*(1+Gn))   ];
elda(1:ne,1:10) = ones(ne,1) * EE1;
femaxi;
%----------------------------------------------------------------------
figure; plot(crd,ur,'-');grid on;xlabel('r [m]');ylabel('u_r [m]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
figure; plot(crdm,Gerr,'-',crdm,Gett,'--',crdm,Gezz,'-.');grid on; 
        xlabel('r [m]');ylabel('\epsilon  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\epsilon_{rr}','\epsilon_{tt}','\epsilon_{zz}');
figure; plot(crdm,Gsrr,'-',crdm,Gstt,'--',crdm,Gszz,'-.');grid on; 
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\sigma_{rr}','\sigma_{tt}','\sigma_{zz}');
figure; plot(crdm,GsTR,'-k',crdm,GsVM,'--k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\sigma_{TR}','\sigma_{VM}');
%======================================================================

%======================================================================
% Open cylinder, internal pressure; plane stress
%----------------------------------------------------------------------
close all; clear all; 
sta = 'pss'; ne = 50; elt = 2;
a = 0.25; b = 0.5; th = 0.5; Gr = 7500; E = 250e9; Gn = 0.33;
Go = 0; Ga = 0; pf=[1 100e6];
Ri = a; Ro = b;
% Element data base 
%     [ E1  E2   E3   Gn12   Gn21   Gn23   Gn32   Gn31  Gn13   G12 ]
EE1 = [ E   E    E    Gn     0      0      Gn     Gn    0      E/(2*(1+Gn))   ];
elda(1:ne,1:10) = ones(ne,1) * EE1;
femaxi;
%----------------------------------------------------------------------
figure; plot(crd,ur,'-');grid on;xlabel('r [m]');ylabel('u_r [m]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
figure; plot(crdm,Gerr,'-',crdm,Gett,'--',crdm,Gezz,'-.');grid on; 
        xlabel('r [m]');ylabel('\epsilon  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\epsilon_{rr}','\epsilon_{tt}','\epsilon_{zz}');
figure; plot(crdm,Gsrr,'-',crdm,Gstt,'--',crdm,Gszz,'-.');grid on; 
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\sigma_{rr}','\sigma_{tt}','\sigma_{zz}');
figure; plot(crdm,GsTR,'-k',crdm,GsVM,'--k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
        legend('\sigma_{TR}','\sigma_{VM}');
%======================================================================

%======================================================================
% Open cylinder, internal pressure, inhomogeneous material
%----------------------------------------------------------------------
close all; clear all; 
sta = 'pss'; ne = 50; elt = 1;
a = 0.25; b = 0.5; E = 250e9; Gn = 0.33; G = E/(2*(1+Gn)); Gr = 7500; 
pf = [1 100e6]; 
Ri = a; Ro = b;
% Element data base 
%     [ E1   E2    E3    Gn12    Gn21    Gn23    Gn32    Gn31   Gn13    G12 ]
EE1 = [ E    E     E     Gn      0       0       Gn      Gn     0       G   ];
EE2 = [ E    10*E  E     Gn/10   0       0       Gn/10   Gn     0       G   ];
elda( 1:25,1:10) = ones(25,1) * EE1; 
elda(26:ne,1:10) = ones(ne-25,1) * EE2;
femaxi;
%----------------------------------------------------------------------
figure; plot(crd,ur,'-k');grid on;xlabel('r [m]');ylabel('u [m]');
        ax=axis;axis([Ri Ro ax(3) ax(4)]);
figure; plot(crdm,Gsrr,'-k',crdm,Gstt,'--k',crdm,Gszz,'-.k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        legend('\sigma_{rr}','\sigma_{tt}','\sigma_{zz}');
        ax=axis;axis([Ri Ro ax(3) ax(4)]);
%======================================================================

%======================================================================
% Disc with central hole, isotropic, plane stress, centrifugal load
%----------------------------------------------------------------------
close all; clear all; 
sta = 'pss'; ne = 50; elt = 2;
a = 0.2; b = 0.5; th = 0.05; Gr = 7500; E = 200e9; Gn = 0.3;Ri = a; Ro = b;
Go = 6*2*pi; Ga = 0;
Ri = a; Ro = b;
% Element data base 
%     [ E1   E2    E3    Gn12    Gn21    Gn23    Gn32    Gn31   Gn13    G12 ]
EE1 = [ E    E     E     Gn      0       0       Gn      Gn     0       E/(2*(1+Gn))   ];
elda(1:ne,1:10) = ones(ne,1) * EE1;
femaxi;
%----------------------------------------------------------------------
figure; plot(crd,ur,'-k');grid on;xlabel('r [m]');ylabel('u [m]');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
figure; plot(crdm,Gsrr,'-k',crdm,Gstt,'--k',crdm,Gszz,'-.k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        legend('\sigma_{rr}','\sigma_{tt}','\sigma_{zz}');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
figure; plot(crdm,GsTR,'-k',crdm,GsVM,'--k');grid on;
        xlabel('r [m]');ylabel('\sigma  [Pa]');
        legend('\sigma_{TR}','\sigma_{VM}');
        ax=axis;axis([0 Ro ax(3) ax(4)]);
%======================================================================

%**********************************************************************
