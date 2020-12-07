This code reproduces the computations described in https://arxiv.org/abs/1810.05885.

The code files that must be executed are the "Step" files. For n=1,...,9, Step{n} produces Res{n}.txt, and requires Res{m}.txt for m<n. .gp files must be executed with Pari/GP, whereas .m files must be executed with Magma.

Additionally, Step5.gp requires the author's package to compute p-adically Galois representations found in the torsion of Jacobians of curves. This package is available at https://github.com/nmascot/LiftTors.

This package is provided in the hope it will be useful, but comes without any guarantee whatsoever.

Nicolas Mascot, December 7, 2020
