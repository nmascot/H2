/* Step 3: further simplify the equation of C */

read("Res2.txt");
print(F); \\ A(y)*x^4+B(y)*x^2+C(y)
A=polcoef(F,4,x);
B=polcoef(F,2,x);
C=polcoef(F,0,x);
D=B^2-4*A*C;
faD=factor(D); \\ D = R*S^2. C projects onto v²=R(u).
printp(faD);
R=faD[1,1];
printp(hyperellratpoints(R,100)); \\ Rat points on v²=R(u)
F=subst(F,y,20+13*y)/13^6;
F=subst(F,x,26/3*x)*3^2/(2*13)^4;
write("Res3.txt",Str("F:=",F,";"));
print(F);
/* (3*y^5 - 6*y^3 + 3*y)*x^4 + (-2*y^8 - 8*y^7 - 4*y^6 + 12*y^5 + 12*y^3 + 4*y^2 - 8*y + 2)*x^2 + (9*y^9 + 36*y^8 - 36*y^7 - 36*y^6 + 18*y^5 + 36*y^4 - 36*y^3 - 36*y^2 + 9*y) */
