%======================================================================
delete('loadincr.m','savedata.m','updaelda.m');

crd0  = [ 0 0; 100 0; 100 100; 0 100 ];

lok = [ ety 1 1 2 3 4 ]; 

elda = [ ety ma A0 0 0 100000 0.3 0 0 0 0 0 0 0 0]; 

pp0   = [ 1 1 0; 1 2 0; 2 1 0; 2 2 0 ];
pp    = [ pp0; 3 2 0; 4 2 0 ];

Gdu = 1; 
pp    = [ pp; 3 1 Gdu; 4 1 Gdu ];

St=0:1:nic+1; Sft=-50:1:nic+1-50;

sada = fopen('savedata.m','w');
s0 = 'ic'; s1 = 'MTp'; s2 = 'Mfi'; s3 = 'Tp'; s4 = 'crd'; s5 = 'eidaB'; 
fprintf(sada,'save(savefile,s0,s1,s2,s3,s4,s5);\n');
fprintf(sada,'SGg(ic) = MTp(3,1)/crd0(3,2); \n');
fprintf(sada,'Sfix(ic) = Mfi(3,1)+Mfi(4,1);\n');
fprintf(sada,'Sfiy(ic) = Mfi(3,2)+Mfi(4,2);\n');
fclose(sada);

mit = 20; ccr = 1e-4;
%======================================================================
