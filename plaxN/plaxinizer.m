%**********************************************************************
%======================================================================
%  Initialisation of nodal values.
%  pe    :  column with prescribed incremental nodal displacements
%  Dp    :  column with iterative nodal displacements
%  Ip    :  column with incremental nodal displacements
%  Tp    :  column with total nodal displacements
%  fe    :  column with external (applied) nodal forces
%  fi    :  column with internal (resulting) nodal forces
%  peC   :  array with prescribed actual displacements
%  feC   :  array with prescribed actual forces
%======================================================================

pe  = zeros(ndof,1); pet = zeros(ndof,1);
p   = zeros(ndof,1); 
Dp  = zeros(ndof,1); Ip  = zeros(ndof,1); Tp  = zeros(ndof,1); 
fe  = zeros(ndof,1); fet = zeros(ndof,1); fi  = zeros(ndof,1); 
fci = zeros(ndof,1); fbi = zeros(ndof,1); fai = zeros(ndof,1); 
peC = zeros(ndof,1); feC = zeros(ndof,1);

%**********************************************************************
