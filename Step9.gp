/* Step 9: Play with the Dokchister resolvents! */

[f,Z,vecs,G,h,C]=read("Res8.txt");

/* Use the resolvents to determine the conj class of rho(Frob_p) */
FindFrob(p)=
{
 my(a,u,ev,M);
 a=Mod(x,f)*Mod(1,p);
 u=a^p*subst(h,'x,a);
 u=trace(u);
 ev=parapply(g->subst(g,'x,u),G);
 ev=select(y->y==0,ev,1);
 if(#ev==0,error("No possible Frobenius"));
 if(#ev>1,print("Found ",#ev," possible Frobenius");return(0));
 ev=ev[1];
 kronecker(p,3)*C[ev]; \\ Twist by -3 so as to get rep attached to autom form / SL(3)
}

/* Example of use:
\r Step9.gp
p=nextprime(10^100)
FindFrob(p)
*/
