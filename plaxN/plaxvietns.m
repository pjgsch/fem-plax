%======================================================================
% Tensile test; linear viscoelastic model; dimensions in mm
%----------------------------------------------------------------------
delete('loadincr.m','savedata.m','updaelda.m');

if ety==10
  Ru = sqrt(A0/pi);
  crd0 = [ 0 0; Ru 0; Ru 100; 0 100 ];
else
  crd0 = [ 0 0; 100 0; 100 100; 0 100 ]; 
end;

lok = [ ety 1 1 2 3 4 ]; 

elda = [ ety ma 0.1 0 0  2.5e5 0.3 12 0 0 0 ];
mm = [ 3.0e6 3.1e-8; 1.4e6 3.0e-7; 3.9e6 3.0e-6;
       5.4e6 2.9e-5; 1.3e6 2.8e-4; 2.3e5 2.7e-3;
       7.6e4 2.6e-2; 3.7e4 2.5e-1; 3.3e4 2.5e+0;
       1.7e4 2.4e+1; 8.0e3 2.3e+2; 1.2e4 2.2e+3 ];

Gdu = 1;
pp0 = [ 1 1 0; 1 2 0; 2 2 0; 4 1 0 ]; 
pp = [ pp0; 3 2 Gdu; 4 2 Gdu ];                

[St,Sft] = mloin(0,0,nic,GDt,'jmp',[0 1]); 

sada = fopen('savedata.m','w');
fprintf(sada,'SGs22(ic+1) = eidaC(1,22);\n');
fclose(sada);

mit = 10; ccr = 0.01; 

%======================================================================
