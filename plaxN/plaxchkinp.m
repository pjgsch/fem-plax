%**********************************************************************
% The input is checked and default parameter values are provided.

prog = 'plax';                                           % program name

% The load incrementation file 'loadincr.m' may be provided elsewhere.
% When it is missing, it is written here.
% The array 'Sft' must be generated with the function 'mloin.m'.

if ~exist('Sft'), Sft(1:1000) = 1; end;

if ~exist('loadincr.m'), 
  loin = fopen('loadincr.m','w');
  fprintf(loin,'if ic==1, peC = pe0*Sft(ic+1);end; \n');
  fprintf(loin,'if ic>1, peC = pe0*(Sft(ic+1)-Sft(ic));end; \n');
  fprintf(loin,'feC = fe0*Sft(ic+1); \n');
  fclose(loin); 
end;

% The data save file 'savedata.m' may be provided elsewhere.
% When it is missing, it is written here.
% When it exists, one line is appended.

if ~exist('savedata.m'), 
  sada=fopen('savedata.m','w');
  fprintf(sada,'Sic(ic+1)=ic;Sit(ic+1)=it;Sti(ic+1)=ti;Snrm(ic+1)=nrm; \n');
  fclose(sada);
else
  sada=fopen('savedata.m','a');
  fprintf(sada,'Sic(ic+1)=ic;Sit(ic+1)=it;Sti(ic+1)=ti;Snrm(ic+1)=nrm; \n');
  fclose(sada);
end;

if ~exist('printout.m'),
prnt=fopen('printout.m','w'); fprintf(prnt,'\n'); fclose(prnt); end;

if ~exist('updaelda.m'),
prnt=fopen('updaelda.m','w'); fprintf(prnt,'\n'); fclose(prnt); end;

% Control parameters are mostly provided individually, but may be
% collectively provided in the array 'cntr'.

if exist('cntr')
  nic = cntr(1); GDt = cntr(2); nl  = cntr(3);
  mit = cntr(4); cnm = cntr(5); ccr = cntr(6);  vrs = cntr(7);
  fat = cntr(8); res = cntr(9); slw = cntr(10); tol = cntr(11);
  deltat = GDt;
end;

if ~exist('matf'), matf  = './mat/defmat'; end;
if ~exist('outf'), outf  = 'defout'; end;

if ~exist('nl'),     nl    = 0;        end;
if ~exist('nic'),    nic   = 1;        end;
if ~exist('mit'),    mit   = 1;        end;
if ~exist('ccr'),    ccr   = 0.01;     end;
if ~exist('deltat'), deltat = 0;       end;
if ~exist('ts'),     ts    = 1;        end;
if ~exist('GDt'),    GDt   = deltat;   end; GDt0 = GDt;
if ~exist('cnm'),    cnm   = 6;        end;
if ~exist('vrs'),    vrs   = 1;        end;
if  exist('fat')~=1, fat   = 0;        end; 
if ~exist('tol'),    tol   = 0;        end;
if ~exist('res'),    res   = 0;        end;
if ~exist('slw'),    slw   = 0;        end;
%if ~exist('eldalivi'), eldalivi = [];  end;
if ~exist('mm'),     mm    = [];       end;
if  exist('tr'),     ntr   = size(tr,1); else, ntr = 0; tr = []; end;

%  Calculation of some parameters from the input data.

nndof = 2;                       % number of nodal degrees of freedom 
nnod  = size(crd0,1);            % number of nodes
ndof  = nnod*nndof;              % number of system degrees of freedom
ne    = size(lok,1);             % number of elements
nenod = size(lok,2)-2;           % number of element nodes 
nedof = nndof*nenod;             % number of element degrees of freedom 
negr  = size(elda,1);            % number of element groups

if nenod==4, neip=4; elseif nenod==8, neip=9; end;

% The loadcase number is set to the initial value 0.

lc = 0;

% The element group number is saved in 'elgr'.

elgr  = lok(:,2);

%**********************************************************************
