/* This steep requires the package available at https://github.com/nmascot/LiftTors to compute Galois representations in the torsion of Jacobians. */
/* This pakage provides in particular install.gp and GalRep.gp */
read("install.gp");
read("GalRep.gp");
read("Res4.txt");

i=Mod('w,'w^2+1);
bar(z)=trace(z+0*i)-z; \\ Complex conjugation
/* Char poly of Frob_p over K(i). Includes twist by -2. */
char1(ap,p)=my(c);c=x^3-ap*x^2+p*bar(ap)*x-p^3;if(kronecker(-2,p)==-1,-subst(c,x,-x),c);
/* Norm of the above from K(i) to K */
char2(ap,p)=my(c);c=char1(ap,p);subst(liftpol(c*bar(c)),'w,0);
p=11;
ap=-7-10*i; \\ From van Geemen & Top
l=3;
e=1024; \\ p-adic accuracy

Chi=char2(ap,p); \\ Charpoly of Frob_11: defines the piece of J[3] we are interested in
if(kronecker(-3,p)==-1,Chi=subst(Chi,x,-x)); \\ Twist by -3

C=[f,g,d0,L,LL,L1,L2,1]; \\ Data describing C
X=GalRep(C,l,p,e,L11,Chi,14); \\ Gal rep found in piece of J[l] where Frob_p acts with charpoly Chi
\\ Working over unram ext of Qp of deg 14 (so as to have l-th roots of 1 needed for pairings)
write("Res5.txt",X);
