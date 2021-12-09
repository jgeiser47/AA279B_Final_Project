Description of each file:

SCRIPTS:
aa279Bnrho.m = This file finds a foundational L2 Halo Orbit.

aa279Bnrho2.m = This file generates the family of Halo Orbits.

aa279Bnrho3.m = This file plots the northern and southern L2 families.

aa279Bnrho4.m = This file does the bulk of the work to find NRHOs from the Halo family.

aa279Bnrho5.m = This file is mainly for making final plots.

FUNCTIONS:
CR3BP_EOM_1.m = CR3BP equations of motion, no state transition matrix, Synodic Frame, Non-Dimensional

CR3BP_EOM_2.m = CR3BP equations of motion, no state transition matrix, Inertial Frame, Dimensional (km,s)

CR3BP_EOM_3.m = includes state transition matrix, Synodic Frame, Non-Dimensional

ZeroCrossingEvent.m = event function that stops Halo Orbit propagation at half-period point when turned on

HowellShooting.m = differential correction method to generate periodic orbits with perpendicular crossings.

CR3BP_Stability.m = Calculates Halo Orbit stability indices. Doesn't work 100% of the time due to unpredictable ordering of eigenvalues.

MAT FILES:
nominal_north_NRHO.mat = Time History for nominal NRHO, non-dimnensional

SMnorthL2.mat = Sun-Mars northern L2 Halo Orbit family. ICs and period for 100 orbits.

SMnorthL2b.mat = Same as SMnorthL2.mat, but used a very slightly different and more precise mass parameter when generating these.

SMnorthL2c.mat = 10 more orbits in the Sun-Mars northern L2 Halo Orbit family. These are very close to Mars.

HaloEigenvalues.mat = 110x6 matrix containing eigenvalues for the monodromy matrix of each Halo Orbit in the family.
