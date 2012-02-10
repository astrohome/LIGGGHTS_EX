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

#ifndef LMP_FIX_SPH
#define LMP_FIX_SPH

#include "fix.h"

namespace LAMMPS_NS {

class FixSPH : public Fix {
 public:
  FixSPH(class LAMMPS *, int, char **);
  ~FixSPH();
  int setmask();
  void init();
  void init_list(int, class NeighList *);
  virtual void post_integrate() {};
  virtual void post_integrate_respa(int, int);

 protected:
  int iarg;
  int kernel_id;
  double h,hinv; //kernel constant
  double **cutsq;
  class NeighList *list;
  int nlevels_respa;
};

}

#endif
