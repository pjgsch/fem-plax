%**********************************************************************
%
% fbibcutback.m
%
% When no convergence is reached within the maximum allowable number
% of iteration steps (mit), the increment can be recycled with 
% reduced prescribed loading.
%
%**********************************************************************

if it==0,
  delete state.mat; 
  delete state0.mat; 
  save state; 
  save('state0','slow');
end;

%----------------------------------------------------------------------
if slw==1 
%----------------------------------------------------------------------
if it>mit
  clear; load state; load state0; slow = 2*slow; save('state0','slow');
else
  if slow>=2, slow = round(slow/2); end; save state; save('state0','slow');
end;
%----------------------------------------------------------------------
end;
%----------------------------------------------------------------------


%**********************************************************************
