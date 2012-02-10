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

FixStyle(wall/sph,FixWallSPH)

#else

#ifndef LMP_FIX_WALL_SPH_H
#define LMP_FIX_WALL_SPH_H

#include "fix_sph.h"

namespace LAMMPS_NS {

class FixWallSPH : public FixSPH {
 public:
  FixWallSPH(class LAMMPS *, int, char **);
  ~FixWallSPH();
  int setmask();
  void init();
  void setup(int vflag);
  void post_force(int vflag);
  void post_force_respa(int vflag, int ilevel, int iloop);

 protected:
  int wallstyle;
  double lo,hi,cylradius;
  double r0,D;
};

}

#endif
#endif
