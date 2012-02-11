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

FixStyle(check/timestep/gran,FixCheckTimestepGran)

#else

#ifndef LMP_FIX_CHECK_TIMESTEP_GRAN_H
#define LMP_FIX_CHECK_TIMESTEP_GRAN_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCheckTimestepGran : public Fix {
 public:
  FixCheckTimestepGran(class LAMMPS *, int, char **);
  int setmask();
  void init();
  void end_of_step();
  double compute_vector(int);

 private:
  class MechParamGran* mpg;
  class PairGranHookeHistory* pg;
  class FixWallGran* fwg;
  class FixPropertyGlobal* Y;
  class FixPropertyGlobal* nu;
  void calc_rayleigh_hertz_estims();
  double rayleigh_time,hertz_time;
  double fraction_rayleigh,fraction_hertz,fraction_skin;
  double fraction_rayleigh_lim,fraction_hertz_lim;
  double vmax; //max relative velocity
  double r_min;
  bool warnflag;
};

}

#endif
#endif
