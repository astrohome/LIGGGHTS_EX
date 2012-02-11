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

FixStyle(insert/rate/region,FixInsertRateRegion)

#else

#ifndef LMP_FIX_INSERT_RATE_REGION_H
#define LMP_FIX_INSERT_RATE_REGION_H

#include "fix_insert_pack.h"

namespace LAMMPS_NS {

class FixInsertRateRegion : public FixInsertPack  {
 public:

  FixInsertRateRegion(class LAMMPS *, int, char **);
  ~FixInsertRateRegion();

 protected:

  virtual void calc_insertion_properties();
  virtual int calc_ninsert_this();
};

}

#endif
#endif
