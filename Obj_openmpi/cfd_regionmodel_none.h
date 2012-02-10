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

#ifndef LMP_CFD_REGIONMODEL_NONE_H
#define LMP_CFD_REGIONMODEL_NONE_H

#include "cfd_regionmodel.h"

namespace LAMMPS_NS {

class CfdRegionmodelNone : public CfdRegionmodel {
 public:
  CfdRegionmodelNone(class LAMMPS *, int, int, char **,class FixCfdCoupling* fc);
  ~CfdRegionmodelNone();
  void init();

 protected:
  class FixPropertyAtom *inRegion;
  class FixPropertyGlobal *outRegion;

  double *inregion;
  double *outregion;
  int nout;

  int nlocal_last;

  virtual void special_settings();
  virtual void rm_update();
};

}

#endif
