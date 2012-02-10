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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(wall/gran/hooke/history,FixWallGranHookeHistory)

#else

#ifndef LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H
#define LMP_FIX_WALL_GRAN_HOOKE_HISTORY_H

#include "fix_wall_gran.h"

namespace LAMMPS_NS {

class FixWallGranHookeHistory : public FixWallGran {
 public:
  FixWallGranHookeHistory(class LAMMPS *, int, char **);
  ~FixWallGranHookeHistory();

 protected:
  virtual void post_create();
  virtual void init_substyle();
  void addHeatFlux(int, double, double area_ratio);
  virtual void compute_force(int ip, double deltan, double rsq,double meff_wall, double dx, double dy, double dz,double *vwall,double *c_history,double area_ratio);
  virtual void addCohesionForce(int &, double &, double &,double area_ratio);
  virtual void deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu);
  virtual void pre_reset_history(int,double*) {}

  int dampflag,cohesionflag,rollingflag;
  double **Yeff,**Geff,**betaeff,**veff,**cohEnergyDens,**coeffRestLog,**coeffFrict,charVel,**coeffRollFrict;

};

}

#endif
#endif
