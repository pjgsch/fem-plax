%**********************************************************************
%
% fbibcol2mat1.m
%
% Transforms columns to matrices for easier interpretation.
%
%**********************************************************************

Mp  = reshape(p,nndof,nnod)';   
MDp = reshape(Dp,nndof,nnod)';   MIp = reshape(Ip,nndof,nnod)';
MTp = reshape(Tp,nndof,nnod)';    
Mfi = reshape(fi,nndof,nnod)';   Mfe = reshape(fe,nndof,nnod)';
Mrs = reshape(rs,nndof,nnod)';

if ntr>=1
MTpT = reshape(TpT,nndof,nnod)';
MfiT = reshape(fiT,nndof,nnod)'; MfeT = reshape(feT,nndof,nnod)';
MrsT = reshape(rsT,nndof,nnod)';
end;

%**********************************************************************
