/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#ifdef SPH_KERNEL_CLASS

    // kernel identifier (a unique integer >= 0)
    // a name for the kernel
    // name of the functions for the kernel, its derivative, and the cutoff are defined
    SPHKernel
    (
        2,
        wendland,
        sph_kernel_wendland,
        sph_kernel_wendland_der,
        sph_kernel_wendland_cut
    )

#else

#ifndef LMP_SPH_KERNEL_WENDLAND
#define LMP_SPH_KERNEL_WENDLAND

namespace SPH_KERNEL_NS {
  inline double sph_kernel_wendland(double s, double h, double hinv);
  inline double sph_kernel_wendland_der(double s, double h, double hinv);
  inline double sph_kernel_wendland_cut();
}

/* ----------------------------------------------------------------------
   Wendland SPH kernel
   h is kernel parameter
   s is distance normalized by h
   0.41778 is 21 over 16pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland(double s, double h, double hinv)
{
    return (0.41778*hinv*hinv*hinv * (1.-0.5*s)*(1.-0.5*s)*(1.-0.5*s)*(1.-0.5*s) * (2.*s+1));
}

/* ----------------------------------------------------------------------
   Derivative of Wendland SPH kernel
   is equal to grad W if multiplied with radial unit vector
   h is kernel parameter
   s is distance normalized by h
   0.41778 is 21 over 16pi
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland_der(double s,double h, double hinv)
{
    return (0.41778*hinv*hinv*hinv * ( (2*(1-0.5*s)*(1-0.5*s)*(1-0.5*s)*(1-0.5*s)) - (2*(1-0.5*s)*(1-0.5*s)*(1-0.5*s)*(2.*s+1)) ));
}

/* ----------------------------------------------------------------------
   Definition of Wendland SPH kernel cutoff in terms of s
   s is normalized distance
------------------------------------------------------------------------- */

inline double SPH_KERNEL_NS::sph_kernel_wendland_cut()
{
    return 2.;
}

#endif
#endif

