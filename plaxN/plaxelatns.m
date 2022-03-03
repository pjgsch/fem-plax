%**********************************************************************
delete('loadincr.m','savedata.m','updaelda.m');

if ety==10
  Ru = sqrt(A0/pi);
  crd0 = [ 0 0; Ru 0; Ru 100; 0 100 ];
else
  crd0 = [ 0 0; 100 0; 100 100; 0 100 ]; 
end;

lok = [ ety 1 1 2 3 4 ]; 

elda = [ ety ma A0 0 0 100000 0.3 0 0 0 0 0 0 0 0]; 

pp0 = [ 1 1 0; 1 2 0; 2 2 0; 4 1 0 ]; 

Gdu = 1; 
pp = [ pp0; 3 2 Gdu; 4 2 Gdu ];  
              
St=0:1:nic; Sft=-50:1:nic-50;

sada = fopen('savedata.m','w');
s0 = 'ic'; s1 = 'MTp'; s2 = 'Mfi'; s3 = 'Tp'; s4 = 'crd'; s5 = 'eidaB'; 
fprintf(sada,'save(savefile,s0,s1,s2,s3,s4,s5);\n');
fprintf(sada,'Sfy(ic+1) = Mfi(4,2)+Mfi(3,2);\n');
fprintf(sada,'SGs22(ic+1) = eidaC(1,22);\n');
fprintf(sada,'SA(ic+1) = crd(3,1)*eidaB(1,3);  \n');
fprintf(sada,'if vrs==3, SA(ic+1) = pi*crd(3,1)^2; end; \n');
fprintf(sada,'SGl(ic+1) = (crd0(3,2)+MTp(3,2))/crd0(3,2); \n');
fclose(sada);

mit = 10; ccr = 1e-4;

%**********************************************************************
