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

#ifdef PAIR_CLASS

PairStyle(born/coul/long,PairBornCoulLong)

#else

#ifndef LMP_PAIR_BORN_COUL_LONG_H
#define LMP_PAIR_BORN_COUL_LONG_H

#include "pair.h"

namespace LAMMPS_NS {

class PairBornCoulLong : public Pair {
 public:
  PairBornCoulLong(class LAMMPS *);
  ~PairBornCoulLong();
  void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart(FILE *);
  void read_restart(FILE *);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);
  void *extract(char *);

 private:
  double cut_lj_global;
  double **cut_lj,**cut_ljsq;
  double cut_coul,cut_coulsq;
  double **a,**rho,**sigma,**c,**d;
  double **rhoinv,**born1,**born2,**born3,**offset;
  double g_ewald;

  void allocate();
};

}

#endif
#endif
