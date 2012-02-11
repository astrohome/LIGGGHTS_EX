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

FixStyle(transportequation/scalar,FixScalarTransportEquation)

#else

#ifndef LMP_FIX_SCALAR_TRANSPORT_EQUATION_H
#define LMP_FIX_SCALAR_TRANSPORT_EQUATION_H

#include "fix.h"

namespace LAMMPS_NS {

class FixScalarTransportEquation : public Fix {
 public:
  FixScalarTransportEquation(class LAMMPS *, int, char **);
  ~FixScalarTransportEquation();
  int setmask();
  void post_create();
  void pre_delete();
  void init();
  void updatePtrs();
  void final_integrate();
  void initial_integrate_respa(int,int,int);
  void initial_integrate(int);
  double compute_scalar();
  bool match_equation_id(const char*);

 private:
  int nlevels_respa;

  char *equation_id;

  class FixPropertyAtom* fix_quantity;
  char *quantity_name;
  class FixPropertyAtom* fix_flux;
  char *flux_name;
  class FixPropertyAtom* fix_source;
  char *source_name;

  //storage capacity - would be thermal capacity for heat conduction
  int capacity_flag;
  class FixPropertyGlobal* fix_capacity;
  double *capacity;
  char *capacity_name;

  double quantity_0;              
  double *quantity;           
  double *flux;       
  double *source;     

};

}

#endif
#endif

