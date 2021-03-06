/* Step 1: get a first model for C for l=3 */

Fa=x*y*(x^2-1)*(y^2-1)*(x^2-y^2+a*x*y); \\ Family of surfaces considered by van Geemen & Top
F=subst(Fa,a,2); \\ Equation of S

F=subst(F,y,l*x)/x^4; \\ Quartic
F=subst(F,x,1/x+1); \\ Cubic
F=subst(F,x,(x/2+l^2)/(1-l^2))*(x+2*l^2)^4/(2*(l+1)*(l-1))^2; \\ Simplify
printp(factor(F)); \\ Equation of generic fibre

m=polcoef(F,3,x); \\ Quadratic-twisting factor l*(l^2-2*l-1)
F=subst(F,x,x/m)*m^2; \\ Untwist to get a monic cubic model

ellDivpol(f,n)= /* f monic cubic -> ell curve E: y²=f(x) -> prod(y-y(P)) for P in E[n] */
{
  my(E,D);
  E=ellinit([0,polcoef(f,2),0,polcoef(f,1),polcoef(f,0)]);
  D=elldivpol(E,n);
  polresultant(D,y^2-f,x);
}

F=ellDivpol(F,3);
F=subst(F,l,x);
write("Res1.txt",Str("F:=",F,";"));
print(F);
/*
-256*x^56 + 6144*x^55 - 62464*x^54 + 333824*x^53 - 859648*x^52 - 120832*x^51 + 7252992*x^50 - 16046080*x^49 - 9891072*x^48 + 90136576*x^47 - 73076736*x^46 - 237805568*x^45 + 420485120*x^44 + 341843968*x^43 - 1165840384*x^42 - 192667648*x^41 + 2178936320*x^40 - 238563328*x^39 - 3063240704*x^38 + 639488000*x^37 + 3412593664*x^36 - 639488000*x^35 - 3063240704*x^34 + 238563328*x^33 + 2178936320*x^32 + 192667648*x^31 - 1165840384*x^30 - 341843968*x^29 + (-288*y^4 + 420485120)*x^28 + (3456*y^4 + 237805568)*x^27 + (-14400*y^4 - 73076736)*x^26 + (14976*y^4 - 90136576)*x^25 + (56160*y^4 - 9891072)*x^24 + (-142848*y^4 + 16046080)*x^23 + (-52992*y^4 + 7252992)*x^22 + (400896*y^4 + 120832)*x^21 + (-55872*y^4 - 859648)*x^20 + (-624384*y^4 - 333824)*x^19 + (134784*y^4 - 62464)*x^18 + (624384*y^4 - 6144)*x^17 + (-55872*y^4 - 256)*x^16 + (16*y^6 - 400896*y^4)*x^15 + (-96*y^6 - 52992*y^4)*x^14 + (-384*y^6 + 142848*y^4)*x^13 + (3232*y^6 + 56160*y^4)*x^12 + (-5424*y^6 - 14976*y^4)*x^11 + (960*y^6 - 14400*y^4)*x^10 - 3456*y^4*x^9 + (960*y^6 - 288*y^4)*x^8 + 5424*y^6*x^7 + 3232*y^6*x^6 + 384*y^6*x^5 - 96*y^6*x^4 - 16*y^6*x^3 + 27*y^8
*/
