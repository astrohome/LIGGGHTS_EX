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

FixStyle(insert/pack,FixInsertPack)

#else

#ifndef LMP_FIX_INSERT_PACK_H
#define LMP_FIX_INSERT_PACK_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertPack : public FixInsert {
 public:

  FixInsertPack(class LAMMPS *, int, char **);
  ~FixInsertPack();
  void init();
  virtual void restart(char *);

 protected:

  virtual void calc_insertion_properties();
  void init_defaults();

  virtual void x_v_omega(int,int&,int&,double&);
  int overlap(int);

  // region to be used for insertion
  class Region *ins_region;
  double region_volume;
  int ntry_mc;

  // the desired volume fraction for packing
  double volumefraction;

};

}

#endif
#endif
