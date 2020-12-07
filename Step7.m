/* Step 7: Decompose PSU(3,3)=SU(3,3) into conjugacy classes */

/* Write a matrix in gp format */
GPmat:=function(M)
	m:=Nrows(M);
	n:=Ncols(M);
	s:="[";
	for i:=1 to m do
		for j:=1 to n do
			s cat:= Sprint(M[i,j]);
			if j lt n then
				s cat:= ",";
			end if;
		end for;
		if i lt m then
			s cat:= ";";
		end if;
	end for;
	s cat:= "]";
	return s;
end function;


F3:=GF(3);
R<x>:=PolynomialRing(F3);
T:=x^2+1;
F9:=ext<GF(3)|T>; /* F9=F3(i) */
A:=IdentityMatrix(F9,3);
G:=SU(A); /* SU(standard Hermitian form) */
C:=ConjugacyClasses(G);
F9<i>:=BaseRing(G); /* Rename generator of F3(i); strangely its name has been changed by the previous line (Magma bug?) */
C:= [[GPmat(g) : g in Conjugates(G,c[3])] : c in C];
Write("Res7.txt","{");
Write("Res7.txt",C);
Write("Res7.txt","}");


