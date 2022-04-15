### Algorithms

Based on the RELIC toolkit we implement:

* About pairing computation
 1. The Elliptic Net algorithm (with lazy reduction).
 2. The improved Elliptic Net algorithm (with lazy reduction) by Chen et al in 2015.
 3. Faster algorithm based on the second algorithm.

We implement the optimal ate pairing on the BLS12-P381 curve and the KSS18-P676 curve. The preset file is in the preset folder, named gmp-pbc-bls381.sh and gmp-pbc-kss676.sh, respectively.

Miller's algorithm is implemented by the authors of this library (Deigo et al), we just run the existed function for comparison.

* About scalar multiplication
 4. Scalar multiplication algorithm with elliptic nets by SubramanyaRao et al in 2019.
 5. Faster algorithm based on the fourth algorithm. 

We implement scalar multiplication on the NIST-P384 curve and the NIST-P521 curve. Similarly, we give their preset files: nistp384-pm.sh and nistp521-pm.sh.

### Build instructions

Instructions for building the library can be found in the [Wiki](https://github.com/relic-toolkit/relic/wiki/Building).

Here we give a direct way to verify our code.

1. Create a target directory : mkdir build
2. Locate the <preset> file in the preset folder and run the following:
cd build
../preset/<preset>.sh ../
make
3. Enter the bin folder and run the corresponding file.

List of our benchmark for implementation of the Ellpitic Net algorithm:

* elgmain_BLS12
* elgmain_KSS18
* scalmul_384
* scalmul_521

Notice that each preset file attached to only one benchmark here.
For example, we want to run the first benchmark, we should follow the three steps and choose the preset file gmp-pbc-bls381.sh.

### Source code

The source code of our algorithms distributed in different folders. Here we give a brief description.

* Pairings: relic-master/src/pp  en_**.c
* Scalar Multiplication: relic-master/src/ep/  scalmul_**.c
* Main function of benchmark: relic-master/bench/  elgmain_**.c scalmul_**.c
* Some other functions: relic-master/src/fpx  tool.c

Notice that the file whose filename contains "tw" represents the computation is on the twisted curves. 