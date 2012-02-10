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

#ifdef BOND_CLASS

BondStyle(fene/expand,BondFENEExpand)

#else

#ifndef LMP_BOND_FENE_EXPAND_H
#define LMP_BOND_FENE_EXPAND_H

#include "stdio.h"
#include "bond.h"

namespace LAMMPS_NS {

class BondFENEExpand : public Bond {
 public:
  BondFENEExpand(class LAMMPS *);
  ~BondFENEExpand();
  void compute(int, int);
  void coeff(int, char **);
  void init_style();
  double equilibrium_distance(int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  double single(int, double, int, int);

 private:
  double TWO_1_3;
  double *k,*r0,*epsilon,*sigma,*shift;

  void allocate();
};

}

#endif
#endif
