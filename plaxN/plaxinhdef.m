%======================================================================
% Example inhomogeneous deformation : deformation tensor
%----------------------------------------------------------------------
delete('loadincr.m','savedata.m');

crd0 = [ 0 0; 1 0; 1 1; 0 1 ]; 
lok = [ ety 1 1 2 3 4 ];

elda = [ ety ma 1 0 0 100e9 0.25 0 0 0 0 0 0 0 0 ];

pp = [ 1 1 0; 1 2 0; 2 2 0; 4 1 0; 4 2 0 ]; 
pp = [ pp; 2 1 0.5; 3 1 0.5; 3 2 0.5];

GDt=1; nic=1; St=[0 1]; Sft=[0 1]; slw=0;
%======================================================================
