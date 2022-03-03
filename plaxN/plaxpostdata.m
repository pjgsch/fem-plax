%**********************************************************************
linelen=70;
line1=[];for i=1:linelen, line1=[line1 '*']; end;
line2=[];for i=1:linelen, line2=[line2 '=']; end;
line3=[];for i=1:linelen, line3=[line3 '-']; end;
%======================================================================
fprintf('%s\n',line1);
fprintf('node   crd0                Dp                  Tp                  crd \n');
for n=1:nnod
  fprintf('%2d    ',n);
  fprintf('%9.2e %9.2e ',crd0(n,1:2));
  fprintf('%9.2e %9.2e ',MDp(n,1:2));
  fprintf('%9.2e %9.2e ',MTp(n,1:2));
  fprintf('%9.2e %9.2e ',crd(n,1:2));
  fprintf('\n');
end;
fprintf('%s\n',line3);
fprintf('node   ext.force           int.force           res.force \n');
for n=1:nnod
  fprintf('%2d    ',n);
  fprintf('%9.2e %9.2e ',Mfe(n,1:2));
  fprintf('%9.2e %9.2e ',Mfi(n,1:2));
  fprintf('%9.2e %9.2e ',Mrs(n,1:2));
  fprintf('\n');
end;

fprintf('%s\n',line2);
gip = 0;

for e=1:ne
  ety = lok(e,1);
  if     ety==3,   neip=4; vrs=2; nenod=4; nedof=8; tol=0;
  elseif ety==11,  neip=4; vrs=1; nenod=4; nedof=8; tol=0;
  elseif ety==10,  neip=4; vrs=3; nenod=4; nedof=8; tol=0;
  elseif ety==100, neip=2; vrs=0; nenod=4; nedof=8; tol=0;
  end;

  fprintf('            11          22          33          12          21          det  \n');

  for ip=1:neip
    gip = gip + 1;

    fprintf('element    %2d    point          %2d  ',e,ip);
    fprintf('ip.crds.   %11.4e %11.4e    \n',eidaB(gip,28),eidaB(gip,29));

    fprintf('deften     ');
    fprintf('%11.4e %11.4e %11.4e %11.4e %11.4e ',eidaB(gip,95:99));
    fprintf('%11.4e \n',eidaB(gip,100));
    fprintf('strain     ');
    fprintf('%11.4e %11.4e %11.4e %11.4e \n',eidaB(gip,17:20));
    fprintf('stress     ');
    fprintf('%11.4e %11.4e %11.4e %11.4e \n',eidaB(gip,21:24));
    fprintf('\n');
  end;

end;


fprintf('%s\n',line1);

%**********************************************************************

