/* Step 8: compute Dokchitser reolvents */

p=11;
h='x^2;

i=Mod(Mod('i,'i^2+1),3);
C=Mod(read("Res7.txt"),3); \\ Conjugacy classes of PSU(3,3)
[f,Zf,vecs]=read("Res6.txt"); \\ Polynomial, its roots, and their indexation by isotropic vectors
vecs=Mod(vecs,3);
NormaliseVec(v)= \\ Normalise element of P^3 so its first non-zero coordinate is 1
{
  my(n=1);
  while(v[n]==0,n++);
  v/v[n];
}

B=vecmax(apply(abs,polroots(f))); \\ bound on modulus of complex roots of f
B=log(vecsum(apply(abs,Vec(h))))+(1+poldegree(h))*log(B); \\ log bound on |h(z)*z'| for z,z' roots of f
B+= log(2)+log(poldegree(f)); \\ A Dok Res of deg N has log coeffs bounded by N*B
e=vecmax(apply(c->#c,C)); \\ Largest deg of a Dok Res
e=ceil(e*B/log(p))+1; \\ Required p-adic accuracy to accurately identify coeffs of Dok Res
print("Accuracy O(",p,"^",e,")");

export(f,p,e);
Zf=parapply(z->Mod(padicappr(f,liftint(z)+O(p^e))[1],p^e),Zf); \\ Hensel-lift roots to this accuracy

/* Dokchitser resolvent for the conj class C, given roots Z, their indexation vecs, and parameter h(x) */
Dok(C,Z,vecs,h)=
{ \\ prod_{s in C}(x-sum_{f(z)=0} h(z)*s(z) )
	my(nC,d,ZR,g,z,hz,v,gz);
	nC = #C;
	d = #Z;
	ZR=vector(nC);
	for(n=1,nC,
		g=C[n];
		for(m=1,d,
			z=Z[m];
			hz=subst(h,'x,z);
			v=vecs[m];
			v=NormaliseVec(g*v);
			gz=Z[select(w->w==v,vecs,1)[1]];
			ZR[n]+=hz*gz;
		);
		ZR[n]='x-ZR[n];
	);
	centerlift(liftpol(factorback(ZR)));
}

export(NormaliseVec,Dok,Zf,vecs,h);
DokRes=parapply(c->Dok(c,Zf,vecs,h),C); \\ Compute resolvents in parallel

/* To index the Dokchitser resolvents, we pick one "nice" representative of ech conj class */

Mstat(M)= /* Beauty rate of a matrix */
{
	my(a,b,c,d);
	for(j=1,3,
		for(k=1,3,
			if(M[k,j]==0 && k!=j,a++); \\ 0 off the diagonal
			if(M[k,j]==0 && k+j!=4,b++); \\ 0 off the anti-diagonal
			\\ Favour the other coefficients in F3 rather than F3(i), esp. if on top-left of dagonal
			if(M[k,j]==1,c+=1+if(k==j,4-j));
			if(M[k,j]==-1,d+1+if(k==j,4-j));
		)
	);
	[a,b,c,d];
}

Creps=centerlift(liftpol(apply(c->c[vecsort(c,Mstat,5)[1]],C))); \\ Pick nicest representative of each conj class
T=Zf[1].mod; \\ The roots live in Fq=Fp[t]/T
write("Res8.txt",[f,[liftall(Mod(Zf,p)),T,p],centerlift(liftpol(vecs)),DokRes,h,Creps]);
