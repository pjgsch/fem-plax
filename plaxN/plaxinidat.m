%**********************************************************************
% plaxinidat.m : part of plax.m
%
% Initialization of data bases from input data.
% Input data are available in 'elda'.
%
%**********************************************************************

function [eida0,eidaB,eidaC,eismB,eismC,elip,neip] = ...
          plaxeidata(ne,elgr,elda,neip,GDt,eldalivi);

% Input
% ne       :  number of elements
% elgr     :  element group (-> g)
% elda     :  element initial data from input
% neip     :  number of element integration points
% GDt      :  time step
% eldalivi :  linear viscoelastic material data (mm)
% Output
% eida0    :  initial element integration point data 
% eidaB    :  begin increment integration point data
% eidaC    :  current element integration point data
% eismB    :
% eismC    :
% elip     :  integration point numbers in elements
% neip     :  number of element integration points

eida0 = zeros(ne*neip,25);     % initial element integration point data 
eidaB = zeros(ne*neip,100);    % begin increment integration point data
eidaC = zeros(ne*neip,100);    % current element integration point data
elip  = zeros(ne,9);
nemalivi = 0; maxmod = 0;

% In a loop over all elements the input data are inserted in the 
% database array 'elda0'.
% For each material model different data are relevant.
% They are available in the input array 'elda'.

gip = 0;

for e=1:ne
  g   = elgr(e);

% ety    :  element type
% mat    :  material type
% thk    :  element thickness
% gf2    :  geometry parameter
% gf3    :  geometry parameter

  ety = elda(g,1);   mat = elda(g,2);
  thk = elda(g,3);   gf2 = elda(g,4);
  gf3 = elda(g,5);  

%----------------------------------------------------------------------

  if     ety==3,   neip=4; 
  elseif ety==11,  neip=4;
  elseif ety==10,  neip=4;
  elseif ety==100, neip=2;
  end;

  if     mat==1,        % linear elastic isotropic behaviour

    emo = elda(g,6);      %
    nuu = elda(g,7);      %
    alf = elda(g,13);     %

  elseif mat==2,        % linear elastic orthotropic behaviour

    e11 = elda(g,6);   
    n12 = elda(g,7);
    e22 = elda(g,8);   
    g12 = elda(g,9);
    e33 = elda(g,10);  
    n31 = elda(g,11);   
    n32 = elda(g,12);  
    ang = elda(g,13);

  elseif mat==3,        % linear Cauchy - Finger

    emo = elda(g,6);   nuu = elda(g,7);    alf = elda(g,13);
    c0  = emo*nuu/((1-2*nuu)*(1+nuu));     c1  = emo/(1+nuu);

  elseif mat==4,

    emo = elda(g,6);   nuu = elda(g,7);    alf = elda(g,13);
    c0  = emo*nuu/((1-2*nuu)*(1+nuu));     c1  = emo/(1+nuu);

  elseif mat==5,

    deltat = GDt;
    emo = elda(g,6);   nuu = elda(g,7); 
    gam = elda(g,8);   NNN = elda(g,9);
    Gsv = elda(g,10);  hhh = elda(g,11);
    hh2 = elda(g,12);  hh3 = elda(g,13);   
    hh4 = elda(g,14);  hh7 = elda(g,15);

  elseif mat==6,

    R    = 8.314;
    T    = 300;
    G    = elda(g,6);
    Gk   = elda(g,7);
    A0   = elda(g,8);
    GDH  = elda(g,9);
    Gt0  = elda(g,10);
    Gm   = elda(g,11);
    Dinf = elda(g,13);
    h    = elda(g,14);
    H    = elda(g,15);

  elseif mat==8,

    nemalivi = nemalivi + 1;
    maxmod   = size(eldalivi,1);
    einf     = elda(g,6);
    nuu      = elda(g,7);
    nmo      = elda(g,8);

  end;

%----------------------------------------------------------------------

  for ip=1:neip
    gip = gip+1;

    elip(e,ip) = gip;

    eida0(gip,1)  = ety;  eida0(gip,2)  = mat;
    eida0(gip,3)  = thk;  eida0(gip,4)  = gf2;
    eida0(gip,5)  = gf3;    

    if     mat==1,

      eida0(gip,6)  = emo;  eida0(gip,7)  = nuu;  eida0(gip,8)  = alf;

    elseif mat==2,

      eida0(gip,6)  = e11;  eida0(gip,7)  = n12;
      eida0(gip,8)  = e22;  eida0(gip,9)  = g12;
      eida0(gip,10) = e33;  eida0(gip,11) = n31;  eida0(gip,12) = n32;
      eida0(gip,13) = ang;

    elseif mat==3,

      eida0(gip,6)  = c0 ;  eida0(gip,7)  = c1;

    elseif mat==4,

      eida0(gip,6)  = c0 ;  eida0(gip,7)  = c1;

    elseif mat==5

      eida0(gip,6)  = emo;  eida0(gip,7)  = nuu;
      eida0(gip,8)  = gam;  eida0(gip,9)  = NNN;
      eida0(gip,10) = Gsv;  eida0(gip,11) = hhh;
      eida0(gip,12) = hh2;  eida0(gip,13) = hh3;  
      eida0(gip,14) = hh4;  eida0(gip,15) = hh7;
      eidaC(gip,70:72) = [1 1 1];
      eidaC(gip,81)    = Gsv;

    elseif mat==6

      eida0(gip,6)  = G;   
      eida0(gip,7)  = Gk;   eida0(gip,8)  = A0;
      eida0(gip,9)  = GDH;  eida0(gip,10) = Gt0;  
      eida0(gip,11) = Gm;   eida0(gip,12) = 0; 
      eida0(gip,13) = Dinf; eida0(gip,14) = h;
      eida0(gip,15) = H;
      eidaC(gip,60:64) = [1 1 1 0 0];
      eidaC(gip,70:74) = [1 1 1 0 0];
%      eidaC(gip,75)    = 1/( 1 + GDt * G / ( A0*exp(GDH/(R*T)) ) );
      eidaC(gip,75)    = 1;
      eidaC(gip,76)    = 0;
%      eidaC(gip,77)    = 1000;
      eidaC(gip,77)    = A0*exp(GDH/(R*T))*Gt0;
      eidaC(gip,78)    = 0; 

    elseif mat==8

      eida0(gip,6) = einf;
      eida0(gip,7) = nuu;
      eida0(gip,8) = nmo;

    end;
  end;   
end;  

% Storage space for linear viscoelastic modal stresses is initialized

eismC = zeros(nemalivi*neip,maxmod*5);
eismB = eismC;

% The initial element integration point data are copied into
% the database arrays 'eidaC' and 'eidaB'.

eidaC(:,1:15) = eida0(:,1:15);   % This must be changed; not needed!
eidaB = eidaC;

%**********************************************************************
