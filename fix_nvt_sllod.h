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

FixStyle(nvt/sllod,FixNVTSlodd)

#else

#ifndef LMP_FIX_NVT_SLODD_H
#define LMP_FIX_NVT_SLODD_H

#include "fix_nvt.h"

namespace LAMMPS_NS {

class FixNVTSlodd : public FixNVT {
 public:
  FixNVTSlodd(class LAMMPS *, int, char **);
  void init();
  void initial_integrate(int);
  void final_integrate();
  void initial_integrate_respa(int, int, int);

 private:
  int nondeformbias;
};

}

#endif
#endif
