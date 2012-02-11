/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/lj126,FixWallLJ126)

#else

#ifndef LMP_FIX_WALL_LJ126_H
#define LMP_FIX_WALL_LJ126_H

#include "fix_wall.h"

namespace LAMMPS_NS {

class FixWallLJ126 : public FixWall {
 public:
  FixWallLJ126(class LAMMPS *, int, char **);
  void precompute(int);
  void wall_particle(int, double);

 private:
  double coeff1[6],coeff2[6],coeff3[6],coeff4[6],offset[6];
};

#endif
}
#endif
