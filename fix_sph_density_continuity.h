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

#ifdef FIX_CLASS

FixStyle(sph/density/continuity,FixSPHDensityContinuity)

#else

#ifndef LMP_FIX_SPH_DENSITY_CONTINUITY_H
#define LMP_FIX_SPH_DENSITY_CONTINUITY_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixSPHDensityContinuity : public FixSPH {
 public:
  FixSPHDensityContinuity(class LAMMPS *, int, char **);
  ~FixSPHDensityContinuity();
  int setmask();
  void init();
  void post_integrate();
  void post_integrate_respa(int, int);
  void reset_dt();

 private:
  double dtdensity;
  double *step_respa;
};

}

#endif
#endif

