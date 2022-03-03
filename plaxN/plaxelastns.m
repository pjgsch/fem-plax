%======================================================================
% Tensile test; elastic models;                       
% We have to limit the axial deformation to prevent the cross-section
% to become zero.
% Total Lagrange (tol=1) for mat=4 for plane strain.
%----------------------------------------------------------------------
delete('loadincr.m','updaelda.m','savedata.m');
Gdu = 1; neip=4; mit=20; ccr=1e-4;
lok = [ 1 1 1 2 3 4 ]; crd0 = [ 0 0; 100 0; 100 100; 0 100 ];
pp0 = [ 1 1 0; 1 2 0; 2 2 0; 4 1 0 ]; pp = [ pp0; 4 2 Gdu; 3 2 Gdu ]; 
sada = fopen('savedata.m','w');
s0 = 'ic'; s1 = 'MTp'; s2 = 'Mfi'; s3 = 'Tp'; s4 = 'crd'; s5 = 'eidaB'; 
fprintf(sada,'save(savefile,s0,s1,s2,s3,s4,s5);\n');
fprintf(sada,'Sfy(ic+1) = Mfi(4,2)+Mfi(3,2);\n');
fprintf(sada,'SGs22(ic+1) = eidaC(1,22);\n');
fprintf(sada,'SA(ic+1) = crd(3,1)*eidaB(1,3);  \n');
fprintf(sada,'if vrs==3, SA(ic+1) = pi*crd(3,1)^2; end; \n');
fprintf(sada,'SGl(ic+1) = (crd0(3,2)+MTp(3,2))/crd0(3,2); \n');
fclose(sada);
if     ii==1; mat=1; nic=250;	     mty=1; st='ps'; ety=3;
elseif ii==2; mat=3; nic=150;	     mty=1; st='ps'; ety=3;
elseif ii==3; mat=3; nic=175; tol=1; mty=1; st='ps'; ety=3;
elseif ii==4; mat=4; nic=150; tol=0; mty=1; st='ps'; ety=11;
elseif ii==5; mat=4; nic=150; tol=1; mty=2; st='pe'; ety=3; end;
elda = [ ety mat 0.1 0 0 100000 0.3 0 0 0 0 0 0 0 0 ];
St=0:1:nic; Sft=-50:1:nic-50;
%======================================================================
