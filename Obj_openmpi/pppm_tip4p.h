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

#ifdef KSPACE_CLASS

KSpaceStyle(pppm/tip4p,PPPMTIP4P)

#else

#ifndef LMP_PPPM_TIP4P_H
#define LMP_PPPM_TIP4P_H

#include "pppm.h"

namespace LAMMPS_NS {

class PPPMTIP4P : public PPPM {
 public:
  PPPMTIP4P(class LAMMPS *, int, char **);

 private:
  void particle_map();
  void make_rho();
  void fieldforce();

  void find_M(int, int &, int &, double *); 
};

}

#endif
#endif
