%**********************************************************************
% As input the next parameters must be provided.
% ne   :  number of elelements
% elt  :  element type
% Ri   :  inner radius
% Ro   :  outer radius
% elda :  element data
%
% Default values for some parameters are set.
% elt  :  element type
% Gr   :  density
% G    :  shear modulus
% Go   :  radial frequency [rad/s]
% pf   :  prescribed nodal forces
% pu   :  prescribed nodal displacements

if ~exist('elt'), elt = 2; end;
if ~exist('Gr'),  Gr = 1; end;
if ~exist('G'),   G = E/(2*(1+Gn)); end;
if ~exist('Go'),  Go = 0; end;
if ~exist('pf'),  pf = []; end;
if ~exist('pu'),  pu = []; end;

npf = size(pf,1); npu = size(pu,1);

% Next columns are needed later to implement boundary conditions properly.

pa = 1:ne+1; if npu>0, pa(pu(:,1)) = -pa(pu(:,1)); end;

% Nodal coordinates 

crd = Ri:((Ro-Ri)/ne):Ro;

% Stiffness parameters are calculated and stored in element data base.

for e=1:ne
  E1   = elda(e,1); E2   = elda(e,2); E3   = elda(e,3);
  Gn12 = elda(e,4); Gn21 = elda(e,5);
  Gn23 = elda(e,6); Gn32 = elda(e,7);
  Gn31 = elda(e,8); Gn13 = elda(e,9);
  G12  = elda(e,10);
  if Gn12~=0 & Gn21==0, Gn21 = Gn12*(E2/E1); end;
  if Gn21~=0 & Gn12==0, Gn12 = Gn21*(E1/E2); end;
  if Gn23~=0 & Gn32==0, Gn32 = Gn23*(E3/E2); end;
  if Gn32~=0 & Gn23==0, Gn23 = Gn32*(E2/E3); end;
  if Gn31~=0 & Gn13==0, Gn13 = Gn31*(E1/E3); end;
  if Gn13~=0 & Gn31==0, Gn31 = Gn13*(E3/E1); end;
  DD = (1-Gn12*Gn21-Gn23*Gn32-Gn31*Gn13-Gn12*Gn23*Gn31-Gn21*Gn32*Gn13);
  A  = (1-Gn32*Gn23)/(DD)*E1;
  B  = (1-Gn31*Gn13)/(DD)*E2;
  C  = (1-Gn12*Gn21)/(DD)*E3;
  Q  = (Gn31*Gn23+Gn21)/(DD)*E1;
  R  = (Gn21*Gn32+Gn31)/(DD)*E1;
  S  = (Gn12*Gn31+Gn32)/(DD)*E2;
  K  = G12;

  Ae = A; Be = B; Qe = Q;
  As = (A*C - R*R)/C; Bs = (B*C - S*S)/C; Qs = (Q*C - R*S)/C;
  if     sta=='psn'
    Ap = Ae; Bp = Be; Qp = Qe; 
    Zs1 = R; Zs2 = S; Ze1 = 0; Ze2 = 0;
  elseif sta=='pss'
    Ap = As; Bp = Bs; Qp = Qs; 
    Zs1 = 0; Zs2 = 0; Ze1 = -(R/C); Ze2 = -(S/C);
  end;
  elda(e,11:17) = [ Ap Bp Qp Zs1 Zs2 Ze1 Ze2 ];
end;

% Initialization to zero value.

K = zeros(ne+1,ne+1); f = zeros(1,ne+1)'; ur = zeros(1,ne+1);

% Element loop for stiffness matrix.

for e=1:ne
  Ap  = elda(e,11); Bp  = elda(e,12); Qp  = elda(e,13);
  Zs1 = elda(e,14); Zs2 = elda(e,15);
  Ze1 = elda(e,16); Ze2 = elda(e,17);
  R1  = crd(e); R2 = crd(e+1);

  if elt==1

  a    = 1/(R2^2 - R1^2);
  Ke11 = (a^2)*(R1^2)* ...
         (...
           (1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) - ...
           (1/2)*(Ap-2*Qp+Bp)*(R2^4)*(R2^(-2)-R1^(-2)) + ...
           (2)*(Ap-Bp)*(R2^2)*log(R2/R1) ...
         );
  Ke12 = (a^2)*(R1*R2)* ...
         (...
           -(1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) + ...
           (1/2)*(Ap-2*Qp+Bp)*(R1^2)*(R2^2)*(R2^(-2)-R1^(-2)) - ...
           (Ap-Bp)*(R2^2+R1^2)*log(R2/R1) ...
         );
  Ke22 = (a^2)*(R2^2)* ...
         (...
           (1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) - ...
           (1/2)*(Ap-2*Qp+Bp)*(R1^4)*(R2^(-2)-R1^(-2)) + ...
           (2)*(Ap-Bp)*(R1^2)*log(R2/R1) ...
         );
  Ke21 = Ke12;
%  fe1  = Gr*Go^2*a*R1*( -(1/4)*(R2^4-R1^4) + (1/2)*R2^2*(R2^2-R1^2) );
%  fe2  = Gr*Go^2*a*R2*(  (1/4)*(R2^4-R1^4) - (1/2)*R1^2*(R2^2-R1^2) );
  fe1  = Gr*Go^2*(1/4)*R1*(R2^2-R1^2);
  fe2  = Gr*Go^2*(1/4)*R2*(R2^2-R1^2);


%  A11  = a^2*R1^2*( (Ap+2*Qp+Bp) + 2*R2^2*(Ap-Bp) );
%  B11  = a^2*R1^2*R2^4*(Ap+Bp-2*Qp);
%  A12  = -a^2*R1*R2*( (Ap+2*Qp+Bp) + R1^2*(Ap-Bp) + R2^2*(Ap-Bp) );
%  B12  = -a^2*R1^3*R2^3*(Ap+Bp-2*Qp);
%  A21  = -a^2*R1*R2*( (Ap+2*Qp+Bp) + R1^2*(Ap-Bp) + R2^2*(Ap-Bp) );
%  B21  = -a^2*R1^3*R2^3*(Ap+Bp-2*Qp);
%  A22  = a^2*R2^2*( (Ap+2*Qp+Bp) + 2*R1^2*(Ap-Bp) );
%  B22  = a^2*R2^2*R1^4*(Ap+Bp-2*Qp);
%  Ke11 = A11*(1/2)*(R2^2 - R1^2) - B11*(1/2)*(R2^(-2) - R1^(-2));
%  Ke12 = A12*(1/2)*(R2^2 - R1^2) - B12*(1/2)*(R2^(-2) - R1^(-2));
%  Ke21 = A21*(1/2)*(R2^2 - R1^2) - B21*(1/2)*(R2^(-2) - R1^(-2));
%  Ke22 = A22*(1/2)*(R2^2 - R1^2) - B22*(1/2)*(R2^(-2) - R1^(-2));
%  fe1  = Gr*Go^2*(1/(4*a))*R1;
%  fe2  = Gr*Go^2*(1/(4*a))*R2;

  elseif elt==2

  a    = 1/(R2 - R1);
  Ke11 = a^2*( ...
               (1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) - ...
               (2)*(Qp+Bp)*(R2)*(R2-R1) + ...
               (Bp)*(R2^2)*log(R2/R1) ...
             );
  Ke12 = a^2*( ...
               -(1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) + ...
               (Qp+Bp)*(R2+R1)*(R2-R1) - ...
               (Bp)*(R1*R2)*log(R2/R1) ...
             );
  Ke22 = a^2*( ...
               (1/2)*(Ap+2*Qp+Bp)*(R2^2-R1^2) - ...
               (2)*(Qp+Bp)*(R1)*(R2-R1) + ...
               (Bp*R1^2)*log(R2/R1) ...
               );
  Ke21 = Ke12;
  fe1  = Gr*Go^2*a * ( (1/12)*R2^4-(1/3)*R2*R1^3+(1/4)*R1^4 );
  fe2  = Gr*Go^2*a * ( (1/12)*R1^4-(1/3)*R1*R2^3+(1/4)*R2^4 );

  end;

  Ke = [ Ke11 Ke12; Ke21 Ke22 ];
  fe = [ fe1 fe2 ]';
  K(e:e+1,e:e+1) = K(e:e+1,e:e+1) + Ke;
  f(e:e+1)       = f(e:e+1)       + fe;
end;

% Taking prescribed displacements to right hand side.
% Eliminating obsolete rows and columns.

if npf>0
  f(pf(:,1)) = f(pf(:,1)) + pf(:,2).*crd(pf(:,1));
end;
if npu>0
  f = f - K(:,pu(:,1))*pu(:,2);
  K(:,pu(:,1)) = []; K(pu(:,1),:) = []; f(pu(:,1)) = [];
end;

% Inversion leads to solution of unknowns.

urr = inv(K)*f;

% Inserting prescribed displacements.

k=1; for j=1:ne+1, if pa(j)>0, ur(j)=urr(k); k=k+1; end;end;
if npu>0, ur(pu(:,1)) = pu(:,2); end;

% Calculating strains and stresses.

for e=1:ne
  R1 = crd(e); R2 = crd(e+1); Rm = (1/2)*(R1+R2); crdm(e) = Rm;
  u1 = ur(e); u2 = ur(e+1);

  if elt==1

  a = 1/(R2^2 - R1^2);
  Gy1 = a*R1*(-Rm + R2^2/Rm); Gy2 = a*R2*(Rm - R1^2/Rm);
  Gy1r = a*R1*(-1 - R2^2/(Rm^2)); Gy2r = a*R2*(1 + R1^2/(Rm^2));

  elseif elt==2

  a = 1/(R2 - R1);
  Gy1 = a*(R2 - Rm); Gy2 = a*(-R1+Rm);
  Gy1r = -a; Gy2r = a;

  end;

  Gerr(e) = Gy1r*u1 + Gy2r*u2;
  Gett(e) = (1/Rm)*(Gy1*u1 + Gy2*u2);
  Gezz(e) = Ze1 * Gerr(e) + Ze2 * Gett(e);
  Gsrr(e) = Ap*Gerr(e) + Qp*Gett(e);
  Gstt(e) = Qp*Gerr(e) + Bp*Gett(e);
  Gszz(e) = Zs1 * Gerr(e) + Zs2 * Gett(e);
  GsTR(e) = max( [abs(Gsrr(e)-Gstt(e)) ...
                  abs(Gsrr(e)-Gszz(e)) abs(Gstt(e)-Gszz(e))] );
  GsVM(e) = sqrt(0.5*((Gsrr(e)-Gstt(e))^2 + ...
                      (Gstt(e)-Gszz(e))^2 + (Gszz(e)-Gsrr(e))^2));
end;

%**********************************************************************
