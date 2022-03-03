%**********************************************************************
% fbibcvnrm.m
%
% Calculates the convergence norm.
%
%**********************************************************************

function [nrm] = fbibcnvnrm(cnm,pu,ppc,prs,Dp,Ip,rs,fi);

nrms=zeros(10,1);

maxDp = max(abs(Dp));  if maxDp<=1e-10, maxDp=1; end;
nrms(1) = (Dp'*Dp)/maxDp;
nrms(2) = max(abs(Dp));
%if size(pu,2)>0, nrms(3) = sqrt((Dp(pu)' * Dp(pu))/max(abs(Dp(pu)))); end;
if size(pu,2)>0, nrms(4) = max(abs(Dp(pu)))/max(abs(Ip(pu))); end;
%if size(pu,2)>0, nrms(5) = max(abs(rs(prs)))/max(max(abs(fi([ppc']))),1e-10); end;

%nrms(6) = (Dp'*Dp)/max( max(abs(Dp)), 1e-10 ); 
%sqrs = sqrt(rs'*rs); 
sqrs = sqrt(rs'*rs); if sqrs<=1e-10, sqrs=1; end;
nrms(6) = min(sqrt(rs(pu)'*rs(pu))/sqrs,sqrt(rs(pu)'*rs(pu)));

%nrms(6) = abs( sqrt(Dp'*Dp)/sqrt(Ip'*Ip) );

%  maxrs   = max(abs(rs(pu))); if maxrs<=1e-10, maxrs=1; end;
%  maxfi   = max(abs(fi)); if maxfi<=1e-10, maxfi=1; end;
%  nrms(6) = sqrt( (rs(pu)' * rs(pu))/maxfi);
if size(pu,2)>0, nrms(7) = max(abs(rs(pu)))/max(abs(rs)); end;
nrms(8) = abs(rs'*Dp);

nrm = nrms(cnm);

%**********************************************************************
