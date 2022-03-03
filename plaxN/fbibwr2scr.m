%**********************************************************************

%if it==2
%  fprintf('ic = %2d ; slow = %2d ; it = %d ; ',nic+1-ic,slow,it)
%  fprintf('nrm = %9.4g \n',nrm);
%else
%  fprintf('                      it = %d ; nrm = %9.4g \n',it,nrm);
%end;

fprintf(' %4d %4d %4d      %9.4g \n',nic+1-ic,slow,it,nrm);

%**********************************************************************
