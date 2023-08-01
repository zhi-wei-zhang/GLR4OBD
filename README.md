# GLR4OBD
GLR-based design for finding optimal biological dose (OBD)

There are 5 files in this folder:

"glr simulations.R" is the main program used to produce simulation results for the ZLY, AGLL and GLR-based designs.

"glr calculations.R" is a supporting source file containing functions for GLR-related calculations. It is referenced in "glr simulations.R".

"boin-et simulations.R" is the program used to produce simulation results for the BOIN-ET design. It makes use of the "boinet" package from the original authors (Takeda et al., 2018).

"bams simulations.R" is the program used to produce simulation results for the BAMS adaptive design. It is a lightly modified version of the program "BAMS-adaptive.R" from the original authors (Lin et al., 2023).

"AMS-cal2.cpp" is a supporting source for the BAMS design. It is referenced in "bams simulations.R".
