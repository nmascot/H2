/* Step 2: Determine canonical embdding of C, and project onto a plane. The resulting equation is simpler. */
QQ:=Rationals();
R<x,y>:=PolynomialRing(QQ,2);
A2:=AffineSpace(R);
load "Res1.txt";
C:=Curve(A2,F);
W:=CanonicalEmbedding(C);
C:=CanonicalImage(C,W);
C:=AffinePatch(C,1);
C:=Elimination(C,[1,2,3,4]);
C:=Curve(C);
_<x1,x2,x3,x4,x5,x6>:=CoordinateRing(C);
C:=Image(map<C->A2 | [x5,x6]>);
print(Genus(C)); /* 7: this is the same curve as before, as opposed to a quotient */
F:=3^3*13*Equation(C);
Write("Res2.txt","{F=");
Write("Res2.txt",F);
Write("Res2.txt",";}");
print(F);
/* 351*x^4*y^5 - 35100*x^4*y^4 + 1285362*x^4*y^3 - 20961720*x^4*y^2 + 148459311*x^4*y - 374594220*x^4 - 8*x^2*y^8 + 864*x^2*y^7 - 34064*x^2*y^6 + 519536*x^2*y^5 + 110400*x^2*y^4 - 23713936*x^2*y^3 - 
    2298974896*x^2*y^2 + 48122966496*x^2*y - 247835237752*x^2 + 208*y^9 - 26624*y^8 + 1124032*y^7 - 779584*y^6 - 1602165344*y^5 + 64797194176*y^4 - 1296163023168*y^3 + 14839155321536*y^2 - 
    93670614701872*y + 252643199407040 */
