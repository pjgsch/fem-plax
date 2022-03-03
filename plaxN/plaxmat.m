%**********************************************************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The proper function(s) for a certain material (mat) is selected.
% The function(s) calculate the stress and the material stiffness.
% Calculated data are available globally.
% Some data are saved in the database.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if     mat==1                      % linear elastic isotropic behaviour
%----------------------------------------------------------------------
%dfie = dfie0;
%dt = dt0;
du = dfie0*Tpe;
[mmM,ccGs,ccGe,vm] = plaxelas1(ipdaC(6:7),vrs,du);        % plaxelas1.m
ccGs=ccGs';ccGe=ccGe';mmGS = zeros(5);
thk = thk0 + thk0*ccGe(3);
%----------------------------------------------------------------------
elseif mat==2                    % linear elastic orthotropic behaviour
%----------------------------------------------------------------------
%dfie = dfie0;
%dt = dt0;
du = dfie0*Tpe;
[mmM,ml,ssl,ccGs,snl,ccGe] = plaxelas2(ipdaC(6:13),vrs,du);%plaxelas2.m
ccGs=ccGs';ccGe=ccGe';mmGS = zeros(5);
thk = thk0 + thk0*ccGe(3);
%----------------------------------------------------------------------
elseif mat==3                                  % linear Cauchy - Finger
%----------------------------------------------------------------------
[mmM,mmGS,ccGs,mF] = plaxelas3(ipda0(6:7),vrs,mF);        % plaxelas3.m
ccGs=ccGs'; 
if (vrs==1 | vrs==2), thk = thk0 * mF(3,3); end;
snmm = 1/2*(mF*mF' - eye(3));
ccGe = [snmm(1,1) snmm(2,2) snmm(3,3) snmm(1,2) snmm(2,1)];
%----------------------------------------------------------------------
elseif mat==4
%----------------------------------------------------------------------
if     tol==0
[mmM,mmGS,ccGs,mF] = plaxelas41(ipda0(6:7),vrs,mF);      % plaxelas41.m
ccGs=ccGs';
if (vrs==1 | vrs==2), thk = thk0 * mF(3,3); end;
snmm = 1/2*(mF'*mF - eye(3));
ccGe = [snmm(1,1) snmm(2,2) snmm(3,3) snmm(1,2) snmm(2,1)];
elseif tol==1
[mmM,mmGS,ccGs,mF] = plaxelas42(ipda0(6:7),vrs,mF);      % plaxelas42.m
ccGs=ccGs';
if (vrs==1 | vrs==2), thk = thk0 * mF(3,3); end;
snmm = 1/2*(mF'*mF - eye(3));
ccGe = [snmm(1,1) snmm(2,2) snmm(3,3) snmm(1,2) snmm(2,1)];
dfie = dfie0;
dt = (thk0/thk)*dt0;
end
%----------------------------------------------------------------------
elseif mat==5
%----------------------------------------------------------------------
if (it>0),
  [ipdaC] = plaxperzs(mF,mFB,GDt,ipda0,ipdaB,ipdaC);      % plaxperzs.m
end; 
[mmM,mmGS] = plaxperzm(mF,mFB,GDt,ipda0,ipdaB,ipdaC);     % plaxperzm.m
eidaC(gip,:) = ipdaC;
ccGe = zeros(5,1)';
ccGs = ipdaC(60:64);
%----------------------------------------------------------------------
elseif mat==6
%----------------------------------------------------------------------
ipdaB(75) = 1;
if (it>0), 
  [ipdaC] = plaxleonsS(mF,mFB,GDt,ipda0,ipdaB,ipdaC);    % plaxleonsS.m
end;
[mmM,mmGS] = plaxleonm(mF,mFB,GDt,ipda0,ipdaB,ipdaC);     % plaxleonm.m
eidaC(gip,:) = ipdaC;
ccGe = zeros(5,1)';
ccGs = ipdaC(65:69);
%----------------------------------------------------------------------
elseif mat==8
%----------------------------------------------------------------------
eldalivi = mm;
vip = neip*(elmalivi-1)+ip;
ipsmB = eismB(vip,:);

ccLut = dfie0*Tpe;
ccGe = [ccLut(1:3); 0.5*(ccLut(4)+ccLut(5)); 0.5*(ccLut(4)+ccLut(5))]';
ipdaC(17:20) = ccGe(1:4);

if (it>0),                                                % plaxviels.m
  [ipdaC,ipsmC] = plaxviels(GDt,ipda0,ipdaB,ipdaC,ipsmB,eldalivi); 
end;
[mmM,mmGS] = plaxvielm(GDt,ipda0,eldalivi,it);            % plaxvielm.m

eidaC(gip,:) = ipdaC;
if (it>0), eismC(vip,:) = ipsmC; end;
ccGs = [ipdaC(21:24) ipdaC(24)];

%----------------------------------------------------------------------
end;

eidaC(gip,16)    = thk;
eidaC(gip,3)     = thk;
eidaC(gip,17:20) = ccGe(1:4); 
eidaC(gip,21:24) = ccGs(1:4); 
eidaC(gip,25)    = dt;

eidaC(gip,90)    = r0;
eidaC(gip,91)    = r;

eidaC(gip,95)    = mF(1,1);
eidaC(gip,96)    = mF(2,2);
eidaC(gip,97)    = mF(3,3);
eidaC(gip,98)    = mF(1,2);
eidaC(gip,99)    = mF(2,1);
eidaC(gip,100)   = det(mF);


%**********************************************************************
