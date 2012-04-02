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

FixStyle(heat/gran,FixHeatGran)

#else

#ifndef LMP_FIX_HEATGRAN_H
#define LMP_FIX_HEATGRAN_H

#include "fix.h"

#define BOLTS   5.670400*pow((double)10,-8)

namespace LAMMPS_NS {   
    

class FixHeatGran : public Fix {
 public:

  FixHeatGran(class LAMMPS *, int, char **);
  ~FixHeatGran();
  int setmask();
  void post_create();
  void init();
  void updatePtrs();
  void post_force(int);
  double compute_scalar();
  void cpl_evaluate(class ComputePairGranLocal *);
  void register_compute_pair_local(class ComputePairGranLocal *ptr);
  void unregister_compute_pair_local(class ComputePairGranLocal *ptr);

  int cpl_is_null()
  {
      if(cpl) return 0;
      return 1;
  }

 private:
     
     struct param {
  double da;
  double db;
  double D;

  double xda;
  double xdb;
  double yda;
  double ydb;

  double R [2];
  double integral;

};
  
  double procedeCalc(double,double,double);
  bool scalar_test(struct param * );
  void calculate_vector(struct param * , int , int );
  double integrate(struct param * );
  double sqr_vect_norm(double , double );
  double third_D_coef(double , double , int );
  double dist(double , double , double , double , double , double );
  double powerCalc(double , double );

  template <int> void post_force_eval(int,int);

  class FixPropertyAtom* fix_temp;
  class FixPropertyAtom* fix_heatFlux;
  class FixPropertyAtom* fix_heatSource;
  class FixPropertyGlobal* fix_conductivity;
  class FixScalarTransportEquation *fix_ste;

  class ComputePairGranLocal *cpl;

  double T0;              
  double *Temp;           
  double *heatFlux;       
  double *heatSource;     
  double *conductivity;
  double nb_int;     //Integration precision
  int cutoff;
  int rad;

  // for heat transfer area correction
  int area_correction;
  double const* const* deltan_ratio;

  class PairGran *pair_gran;
  int history_flag;
};

}

#endif
#endif