%======================================================================
% newtonS : Newton stress iteration with subincrements

resGl = 1 - Gl*(Gdt1*GG + 1);
resD  = DJ - D + Gdt2*DD;
nrmGl = abs(resGl);
nrmD  = abs(resD);

nit = 1; mnit = 50;
epsGl = 1e-4; epsD = 1e-4;
nrmGl = 1e3;  nrmD = 1e3; 

while ((nrmGl>epsGl | nrmD>epsD) & nit<=mnit)
%while (nit<=mnit)

  cdd1 = DD/(sqrt(3)*Gt0);
  cdd3 = DD - h*GsVM/(sqrt(6)*Dinf*Gh);

  if (GsVM>0)
  ce11a = exp(GsVM/(sqrt(3)*Gt0));
  ce12a = exp(-GsVM/(sqrt(3)*Gt0));
  ce1a  = Gh + GsVM*Gh*(ce11a+ce12a)/(ce11a-ce12a);
  else
  ce1a = 0;
  end;

  dfdGl = Gdt1*GG + 1 - Gl*Gdt1*GG*(1 - GsVM/(sqrt(3)*Gt0));
%  dfdGl = Gdt1*GG + 1 - Gl*Gdt1*GG/(Gh)*ce1a;
  dfdD  = Gl*Gdt1*GG;
  dgdGl = - Gdt2*cdd1*GsVM;
  dgdD  = 1 - Gdt2*cdd3;

  mat  = [dfdGl dfdD ; dgdGl dgdD];
  res  = [resGl resD]';
  sol  = mat\res;
  dsol = mat\(mat*sol-res);
  soll = sol-dsol;

  dGl = soll(1);    dD = soll(2);
  Gl  = Gl + dGl;    D = D + dD; 

  mCpn    = (1-Gl)*mtCn + Gl*mCpnJ;
  mbtBen  = mtUn*inv(mCpn)*mtUn';
%  mbtBen  = mbtBen ./ (det(mbtBen)^(1/3));
  mbtBend = mbtBen - 1/3*trace(mbtBen)*mI;

  mbsd  = G*mbtBend;

  GsVM  = sqrt(3/2*trace(mbsd*mbsd));

  A     = A0*exp(GDH/(R*T) + Gm*p/Gt0 - D);
  if GsVM>1.0e-10, Gh = A*GsVM/(sqrt(3)*sinh(GsVM/(sqrt(3)*Gt0))); end;
  GG    = G/Gh;
  DD    = h*(1-D/Dinf)*GsVM/(sqrt(6)*Gh);

  resGl = 1 - Gl*(Gdt1*GG + 1);
  resD  = DJ - D + Gdt2*DD;

  nrmGl = abs(dGl);
  nrmD  = abs(dD);
  nit = nit + 1;

end;

%======================================================================
