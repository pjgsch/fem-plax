%**********************************************************************
% This function writes the file 'savedata.m' which is called at the end
% of an increment to save data for postprocessing.
% The input is a combination of a two strings or a string and a matrix.
% Examples :
% msada('s','crd','s','MTp');     saves 'crd' and 'MTp' to mat/.
% msada('dsp',[3 1; 4 1]);        displacements dof 1 in node 3 and 4
% msada('rfc',[3 1; 4 1]);        reaction forces dof 1 in node 3 and 4
% msada('sts',[1],'stn',[1]);     stress and strain in ip 1
% msada('eld',[2 6]);             element 2, parameter 6
% msada('ipd',[2 1 20; 2 1 22]);  element 2, ip 1, parameter 22
%**********************************************************************

function msada(s1,as1,s2,as2,s3,as3,s4,as4,s5,as5,s6,as6);

arg=0;

sada = fopen('savedata.m','a');

while arg<nargin

arg=arg+2;

if     arg==2,  s=s1;  as=as1;
elseif arg==4,  s=s2;  as=as2;
elseif arg==6,  s=s3;  as=as3;
elseif arg==8,  s=s4;  as=as4;
elseif arg==10, s=s5;  as=as5;
elseif arg==12, s=s6;  as=as6;
elseif arg==14, s=s7;  as=as7;
elseif arg==16, s=s8;  as=as8;
elseif arg==18, s=s9;  as=as9;
elseif arg==20, s=s10; as=as10;
end;

if s=='s'
%  fprintf(sada,'save(savefile,''-append'',''%s'');\n',as);
  fprintf(sada,'save(savefile,''crd'');\n');
%  fprintf(sada,'save eval(savefile) crd eldaC \n');
end;

if s=='dsp' | s=='STp'
  sas=size(as,1);
  for i=1:sas
    n=as(i,1); d=as(i,2);
    sn=num2str(n); sd=num2str(d);
    fprintf(sada,'Su');
    fprintf(sada,sn);
    fprintf(sada,sd);
    fprintf(sada,'(ic+1) = MTp(%d,%d); \n',n,d);
    fprintf(sada,'STp');
    fprintf(sada,sn);
    fprintf(sada,sd);
    fprintf(sada,'(ic+1) = MTp(%d,%d); \n',n,d);
  end;
end;

if s=='rfc' | s=='Sfi'
  sas=size(as,1);
  for i=1:sas
    n=as(i,1); d=as(i,2);
    sn=num2str(n); sd=num2str(d);
    fprintf(sada,'Sfi');
    fprintf(sada,sn);
    fprintf(sada,sd);
    fprintf(sada,'(ic+1) = Mfi(%d,%d); \n',n,d);
  end;
end;

if s=='efc'
  sas=size(as,1);
  for i=1:sas
    n=as(i,1); d=as(i,2);
    sn=num2str(n); sd=num2str(d);
    fprintf(sada,'Sfe');
    fprintf(sada,sn);
    fprintf(sada,sd);
    fprintf(sada,'(ic+1) = Mfe(%d,%d); \n',n,d);
  end;
end;

if s=='stn'
  sas=size(as,1);
  for i=1:sas
    e=as(i,1); 
    se=num2str(e);
    fprintf(sada,'SGe');
    fprintf(sada,se);
    fprintf(sada,'(ic+1) = eldaB(%d,6); \n',e);
  end;
end;

if s=='sts'
  sas=size(as,1);
  for i=1:sas
    e=as(i,1); 
    se=num2str(e);
    fprintf(sada,'SGs');
    fprintf(sada,se);
    fprintf(sada,'(ic+1) = eldaB(%d,7); \n',e);
  end;
end;

if s=='ipd'
  sas=size(as,1);
  for i=1:sas
    el=as(i,1);ip=as(i,2);id=as(i,3);
    sel=num2str(el);sip=num2str(ip);sid=num2str(id);
    fprintf(sada,'S');
    fprintf(sada,sel);
    fprintf(sada,sip);
    fprintf(sada,sid);
    fprintf(sada,'(ic+1) = eidaB(elip(%d,%d),%d); \n',el,ip,id);
  end;
end;

if s=='eld'
  sas=size(as,1);
  for i=1:sas
    el=as(i,1);id=as(i,2);
    sel=num2str(el);sid=num2str(id);
    fprintf(sada,'S');
    fprintf(sada,sel);
    fprintf(sada,sid);
    fprintf(sada,'(ic+1) = eldaB(%d,%d); \n',el,id);
  end;
end;



end;

fclose(sada);

%**********************************************************************

    
  
