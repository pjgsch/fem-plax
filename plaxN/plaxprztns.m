%======================================================================
delete('loadincr.m','savedata.m');

if     ety==3,  crd0  = [0 0; 100 0; 100 100; 0 100 ];
elseif ety==11, crd0  = [0 0; 100 0; 100 100; 0 100 ];
elseif ety==10, crd0  = [0 0; sqrt(A0/pi) 0; sqrt(A0/pi) 100; 0 100 ]; 
end;

lok = [ ety 1 1 2 3 4 ]; 

elda = [ ety ma A0 0 0 1800 0.372 0.001 3 37 -200 500 700 800 30000 ]; 

pp0 = [ 1 1 0; 1 2 0; 2 2 0; 4 1 0 ];

Gdu = 1;
pp = [ pp0; 4 2 Gdu; 3 2 Gdu ]; 

[St,Sft,nic,GDt,tend] = mloin(0,0,nic,GDt,'pol',[0 0 nic*GDt 40]); 
%[St,Sft,nic,GDt,tend] = mloin(0,0,nic,GDt,['400/(400*' num2str(GDt) ')*t'],[0 0]); 

msada('s','crd','s','MTp'); 
msada('dsp',[3 2],'rfc',[3 2; 4 2],'ipd',[1 1 61]);

mit = 10; ccr = 1e-5; 

%======================================================================
