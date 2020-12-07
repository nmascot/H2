/* Step 6: Find the hermitian form fixed by Galois, and simplify the polynomial corresponding to isotropic vectors */

\\ The first 3 functions establish a correspondence Fl^d <-> {0,1,...,l^d-1}
\\ and give the action of a matrix under this correspondence

c2i(c,l)= \\ Conversion coordinates -> index in an Fl-space of finite dim
{
  my(i);
  i=subst(Pol(c),'x,l);
  if(i,i,l^#c);
}

i2c(i,l,d)= \\ Conversion index -> corrdinates in an Fl-space of dim d
{
  my(j=i,c=vector(d)~);
  for(n=1,d,
    c[d+1-n] = j%l;
    j = j\l
  );
  c;
}

ActOni(M,i,l)= \\ Action of the matrix M on Fl^d, in terms of indices
{
  my(c,d);
  d = #M;
  c = i2c(i,l,d);
  c = lift(Mod(M*c,l));
  c2i(c,l);
}

\\ Normalise an element of proj space so its first non-zero coordinate is 1
NormaliseVec(v)=
{
	my(n=1);
	while(v[n]==0,n++);
	v/v[n];
}

X=read("Res5.txt");
l=3;
p=11;
d=6;
Chi=sum(i=0,d,x^i); \\ Charpoly of Frob_p
MFrob=matcompanion(Chi);
N=(l^d-1)/8; \\ 1/8 * #F3[Frob_p]*
until(a^4==-1,a=Mod(Mod(random(2*x^d),Chi),l)^N); \\ Get 8th root of 1 in F3[Frob].
M9=centerlift(subst(liftpol(a),x,MFrob)); \\ This has to generate F3(i)*.

/* Gather the roots along F3(i)*-orbits, i.e. projectivise */
Z=vector(N,i,1);
tags=vector(N);
done=vector(l^d-1);
iZ=1;
{
for(i=1,l^d-1,
	if(done[i],next);
	tags[iZ]=j=i;
	while(!done[j],
		Z[iZ] *= X[2][j];
		done[j]=1;
		j = ActOni(M9,j,l);
	);
	iZ++;
);
}
/* Determine polynomial describing the projective representation */
F=factorback(apply(z->'x-z,Z));
F=bestappr(liftpol(F));
F=subst(F,variables(F)[2],0);

/* It has 2 factors, one for isotropc vectors, and the other for non-isotropic vectors */
[F1,F2]=factor(F)~[1,];
/* Discard non-isotropic vectors and the corresponding roots */
Z=Mod(Z,p); \\ Work mod p as p-adic accuracy not requred anymore
Isot=select(z->subst(F1,x,z)==0,Z,1);
Z1=vecextract(Z,Isot);
tags=vecextract(tags,Isot);

/* Find F3(i)-basis of the space F3^6, with explicit change of coordinates */
MI=Mod(M9^2,l); \\ Matrix of i acting of F3^6
B=matrix(d,d); \\ Get an F3-basis of F3^6 of the form v1,i*v1,v2,i*v2,v3,i*v3
{ \\ -> v1,v2,v3 is an F3(i)-basis
until(matdet(B),
	for(j=1,d/2,
		B[,2*j-1]=vector(d,i,random(Mod(1,l)))~;
		B[,2*j]=MI*B[,2*j-1];
	);
);
}
i=ffgen(Mod('i^2+1,3));
bar(z)=trace(z)-z; \\ Conjugation on F3(i)
C=Mat([1,i]);
C=matconcat(matdiagonal(vector(d/2,j,C)));
C=C*B^-1; \\ Transition matrix from canonical F3-basis of F3^6 to the F3(i)-basis v1,v2,v3
Isot=apply(n->C*i2c(n,l,d),tags); \\ Apply this change of coords to express the isotropic vectors on this F3(i)-basis

/* Find unique Hermitian form (up to scaling) for which these vectors are isotropic */
\\ By using indeterminate coefficients: A*N(x1)+...+B*Tr(x1*x2bar)+...
K=matrix(28,9);
{for(k=1,28,
	v=Isot[k];
	K[k,1]=norm(v[1]);
	K[k,2]=norm(v[2]);
	K[k,3]=norm(v[3]);
	K[k,4]=trace(v[1]*bar(v[2]));
	K[k,5]=trace(i*v[1]*bar(v[2]));
	K[k,6]=trace(v[1]*bar(v[3]));
	K[k,7]=trace(i*v[1]*bar(v[3]));
	K[k,8]=trace(v[2]*bar(v[3]));
	K[k,9]=trace(i*v[2]*bar(v[3]));
);}
K=matker(K);
print("Dim of space of Hermitian forms having these isotropic vectors: ",#K);
K=K[,1]; \\ Coeffs of the Hermitian form
A=matrix(3,3); \\ Matrix of the Hermitian form w.r.t v1,v2,v3
for(j=1,3,A[j,j]=K[j]);
A[1,2]=K[4]+i*K[5];
A[1,3]=K[6]+i*K[7];
A[2,3]=K[8]+i*K[9];
A[2,1]=bar(A[1,2]);
A[3,1]=bar(A[1,3]);
A[3,2]=bar(A[2,3]);

/* Orthonormalise basis by Gram-Schmidt */
dotp(v,w)=v~*A*apply(bar,w);
B=matid(3);
OB=vector(3);
{
for(r=1,3,
	until(Nv, \\ Find non-isotropic vector in span(B)
		v=B*vector(#B,j,random(i))~;
		Nv=dotp(v,v);
	);
	until(norm(s)*Nv==1,s=random(i)); \\ Rescale so it has norm 1
	v*=s;
	OB[r]=v; \\ Append to new basis
	B=matimage(matconcat(apply(w->w-dotp(w,v)*v,Vec(B)))); \\ Restrict to orthogonal
);
}
B=matconcat(OB)^-1; \\ Reverse change of coordinates

/* Express the isotropic vectors on the orthonormal basis -> indexation of roots by isotropic elements of proj space w.r.t. standard Hermitian form */
tags=apply(v->NormaliseVec(B*v),Isot);

/* Simplify polynomial */
[f1,g]=polredabs(F1,1); \\ f1 = simpler polynomial, g = change of variable
g=liftpol(modreverse(g));
z1=apply(z->subst(g,x,z),Z1); \\ Apply change of variable to roots -> roots of simpler polynomial, indexed as above
write("Res6.txt",[f1,z1,tags]);
\\nfisincl(F2,polcompositum(f1,f1)[2]) \\ Check that F1 and F2 have the same splitting field
