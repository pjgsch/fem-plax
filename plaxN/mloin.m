%**********************************************************************
% [St,Sft,nic,tinc,tend] = mloin(tstart,tend,nic,tinc,func,arr)
%
% St      =  (pseudo) time
% Sft     =  function value
% nic     =  number of increments
% tinc    =  time increment
% tstart  =  start time
% tend    =  end time
% func    =  excitation function
% arr     =  array with data points, coefficients or offset
%
% Input parameters :
%  tend < tstart   ->  tend = nic * tinc
%  nic = 0         ->  nic  = round( (tend-tstart)/tinc ) + 1
%  tinc = 0        ->  tinc = (tend-tstart)/nic
%
% Examples :
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'saw',[20 0 0 0.1]); plot(St,Sft);grid on;
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'blk',[2 1 -1 0 0]); plot(St,Sft);grid on;
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'jmp',[0 1 4 2]); plot(St,Sft);grid on;
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'pol',[0 0 2 10 3 10 4 0]); plot(St,Sft);grid on;
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'sin(2*t)'); plot(St,Sft);grid on;
%  [St,Sft,nic,tinc,tend] = mloin(0,10,100,0,'sin(2*t)',[0 2]); plot(St,Sft);grid on;
%
%**********************************************************************


function [Sti,Sft,nic,tinc,tend] = mloin(tstart,tend,nic,tinc,func,arr);

if nargin<6, arr = [0 0]; end;

if tend<=tstart, tend = nic*tinc; end;
if nic==0,       nic = round( (tend-tstart)/tinc ) + 1; end;
if tinc==0,      tinc=(tend-tstart)/nic; end;

i = 0; ft = 0;

if strcmp(func,'saw');

aa = arr(1); hshift = arr(2); vshift = arr(3); coef = arr(4); 
bb = 2*aa; cc = 0;
for t=tstart:tinc:tend
  i=i+1; t = t+hshift; Sti(i) = t; 
  if i==1, ccc = coef; 
  elseif rem(i-1,aa)==cc, ccc = -ccc; cc = arr(1); aa = bb; end;
  ft = ft + ccc*tinc;
%  if     t>0     & t<=aa,          ft = ft + coef*tinc; 
%  elseif t>aa    & t<=3*aa,   ft = ft - coef*tinc;  
%  elseif t>3*aa  & t<=5*aa,   ft = ft + coef*tinc; 
%  elseif t>5*aa  & t<=7*aa,   ft = ft - coef*tinc; 
%  elseif t>7*aa  & t<=9*aa,   ft = ft + coef*tinc; 
%  elseif t>9*aa  & t<=11*aa,  ft = ft - coef*tinc; 
%  elseif t>11*aa & t<=13*aa,  ft = ft + coef*tinc; 
%  end; 
  Sft(i) = ft + vshift;
end;

elseif strcmp(func,'blk')

hshift = arr(4); vshift = arr(5);
for t=tstart:tinc:tend
  tt = t + hshift;
  i=i+1; Sti(i) = tt; aa=arr(1); 
  amp1 = arr(2) + vshift; amp2 = arr(3) + vshift;
  if     tt>=0     & tt<2*aa,   ft = amp1; 
  elseif tt>=2*aa  & tt<4*aa,   ft = amp2;  
  elseif tt>=4*aa  & tt<6*aa,   ft = amp1; 
  elseif tt>=6*aa  & tt<8*aa,   ft = amp2; 
  elseif tt>=8*aa  & tt<10*aa,  ft = amp1; 
  elseif tt>=10*aa & tt<12*aa,  ft = amp2; 
  elseif tt>=12*aa & tt<14*aa,  ft = amp1; 
  end; 
  Sft(i) = ft;
end;

elseif strcmp(func,'jmp'), 

njmp = round(size(arr,2)/2); j = 1; ft = 0; amp = 0;
for t=tstart:tinc:tend
  i = i+1; Sti(i) = t; 
  if j<=njmp & t>=arr(2*j-1), amp = arr(2*j); j = j+1; end; 
  ft = amp;
  Sft(i) = ft;
end;

elseif strcmp(func,'pol')

npts=round(size(arr,2)/2); j=1; ft=0;
while (j+1)<=npts
  t1 = arr(2*j-1); t2 = arr(2*(j+1)-1);
  f1 = arr(2*j);   f2 = arr(2*(j+1));
  if t2==t1
    t = t2; ft = f2; i=i+1; Sti(i)=t; Sft(i)=ft;
  else
    for t=t1:tinc:(t2)
      i=i+1; Sti(i)=t;
      ft = f1+((f2-f1)/(t2-t1))*(t-t1);
      Sft(i) = ft;
    end;
  end;
  j=j+1;
end;


else

hshift=arr(1); vshift=arr(2);
for t=tstart:tinc:tend
  i = i+1;
  t = t+hshift;
  Sti(i) = t;
  ft = vshift + eval(func);
  Sft(i) = ft;
end;
end;

nn = size(Sti,2);
Sdft = diff(Sft);

%figure; plot([0 Sti],[0 Sft]); grid on;
%figure; plot(Sti(1:nn-1),Sdft); grid on;

%**********************************************************************
