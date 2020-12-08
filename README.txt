This code reproduces the computations described in https://arxiv.org/abs/1810.05885.

The code files that must be executed are the "Step" files. For n=1,...,9, Step{n} produces Res{n}.txt, and requires Res{m}.txt for m<n. .gp files must be executed with Pari/GP, whereas .m files must be executed with Magma.

Additionally, Step5.gp requires the author's package to compute p-adically Galois representations found in the torsion of Jacobians of curves. This package is available at https://github.com/nmascot/LiftTors.

The file "Final results.txt" contains the data obtained at Step 8, as well as some extra data mentioned in the article.

This code and these data are provided in the hope that they will be useful, but come without any guarantee whatsoever.

Nicolas Mascot, December 8, 2020
