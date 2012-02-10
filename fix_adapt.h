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

#ifdef FIX_CLASS

FixStyle(adapt,FixAdapt)

#else

#ifndef LMP_FIX_ADAPT_H
#define LMP_FIX_ADAPT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixAdapt : public Fix {
 public:
  FixAdapt(class LAMMPS *, int, char **);
  ~FixAdapt();
  int setmask();
  void post_create();
  void pre_delete();
  void init();
  void pre_force(int);

 private:
  int nadapt;
  int *which;
  char **pair,**param,**var;
  int *ilo,*ihi,*jlo,*jhi;

  int *ivar;
  class Pair **pairptr;
  int *pairindex;
  int *awhich;

  class FixPropertyAtom *fppat;
  char fixid[100];
  int atomstyle_flag;
};

}

#endif
#endif
