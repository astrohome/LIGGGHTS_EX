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

FixStyle(store/coord,FixStoreCoord)

#else

#ifndef LMP_FIX_STORE_COORD_H
#define LMP_FIX_STORE_COORD_H

#include "fix.h"

namespace LAMMPS_NS {

class FixStoreCoord : public Fix {
 public:
  FixStoreCoord(class LAMMPS *, int, char **);
  ~FixStoreCoord();
  int setmask();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();

 private:
  double **xoriginal;         // original coords of atoms
};

}

#endif
#endif
