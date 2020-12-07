/* Step 4: determine some Riemann-Roch spaces on C */

QQ:=Rationals();
R<x,y>:=PolynomialRing(QQ,2);
load "Res3.txt";

RRbasis:=function(D) /* D -> basis of L(D) */
  L,l:=RiemannRochSpace(D);
  L:=[l(Basis(L)[i]) : i in [1..Dimension(L)]];
  return L;
end function;

RRGP:=procedure(D0,E1,E2,path) /* Given divisors D0,E1,E2 on C, write data needed at next step in gp format */
	f:=Equation(AffinePatch(Curve(D0),1));
  _<x,y>:=Parent(f);
  Write(path,"{f=");
  Write(path,f);
  Write(path,";g=");
  Write(path,Genus(Curve(D0)));
  Write(path,";d0=");
  Write(path,Degree(D0));
  Write(path,";L=");
	L:=RRbasis(D0);
	_<x,y>:=Parent(L[1]);
  Write(path,L);
  Write(path,";LL=");
  Write(path,RRbasis(2*D0));
  Write(path,";L1=");
  Write(path,RRbasis(2*D0-E1));
  Write(path,";L2=");
  Write(path,RRbasis(2*D0-E2));
  Write(path,";}");
end procedure;

Lp:=function(f,p) /* Local L factor at p */
  return LPolynomial(Curve(AffinePlane(GF(p)),f));
end function;

LpGP:=procedure(f,p,path) /* Write local L factor at p in gp format */
  P:=Lp(f,p);
  _<x>:=Parent(P);
  P:=Parent(P)!Reverse(Coefficients(P));
  Write(path,"{L");
  Write(path,p);
  Write(path,"=");
  Write(path,P);
  Write(path,";}");
end procedure;

A2:=AffineSpace(R);
C:=Curve(A2,F);
G:=AutomorphismGroup(C); /* Check shape of phi_C */
print(Order(G));
print(G.2);
KC<x,y>:=FunctionField(C);
P:=Support(Divisor(x))[1];
Q:=Support(Divisor(y))[2];

RRGP(9*P+7*Q,6*P+3*Q,5*P+4*Q,"Res4.txt");
LpGP(F,11,"Res4.txt");
