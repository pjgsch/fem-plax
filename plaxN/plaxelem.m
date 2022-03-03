%**********************************************************************

  if     ety==3,   neip=4; vrs=2; nenod=4; nedof=8; %tol=0;
  elseif ety==11,  neip=4; vrs=1; nenod=4; nedof=8; %tol=0;
  elseif ety==10,  neip=4; vrs=3; nenod=4; nedof=8; %tol=0;
  elseif ety==100, neip=2; vrs=0; nenod=4; nedof=8; %tol=0;
  end;

  for ip=1:neip
    gip = gip + 1;

    ipda0 = eida0(gip,:); ipdaB = eidaB(gip,:); ipdaC = eidaC(gip,:);

    if ety==100
      [em,ef,Gdn,Gdt,eidaC] = ...                       % plaxczelem1.m
      plaxczelem1(em,ef,eidaC,d0,Gfn,Gft,Tnmax,Ttmax,ec0,ec,ip,gip,Tpe); 
    else
      [thk0,thkB,thk,dt0,dtB,dt,dfie0,dfieB,dfie,mF,mFB,r0,r] = ...
      plaxgeom(ec0,ecB,ec,psi,psidksi,vrs,ipda0,ip,nenod); % plaxgeom.m

      vole0 = vole0 + dt0 * ipwf(ip); voleC = voleC + dt * ipwf(ip);

      plaxmat;                                              % plaxmat.m

      em = em + dfie' * (mmM + mmGS) * dfie * thk * dt * ipwf(ip);
      ef = ef + dfie' * ccGs' * thk * dt * ipwf(ip);
    end;
  end;



%**********************************************************************
