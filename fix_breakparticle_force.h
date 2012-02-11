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

FixStyle(breakparticle/force,FixBreakparticleForce)

#else

#ifndef LMP_FIX_BREAKPARTICLE_FORCE_H
#define LMP_FIX_BREAKPARTICLE_FORCE_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixBreakparticleForce : public FixInsert {
 public:

  FixBreakparticleForce(class LAMMPS *, int, char **);
  ~FixBreakparticleForce();

  void post_create();
  void pre_delete();

  virtual int setmask();
  virtual void init();
  void init_defaults();
  int calc_ninsert_this();
  virtual void end_of_step();

 protected:

  // inherited functions
  virtual void calc_insertion_properties();
  virtual void pre_insert();
  int overlap(int);
  void x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this);

  // functions declared in this class
  inline void generate_random(double *pos, double rad_broken,double rad_insert);
  void print_stats_breakage_during();

  // per breakage flag
  class FixPropertyAtom *fix_break;

  // template holding data of the fragments
  class FixTemplateMultiplespheres *fix_fragments;

  // force treshold and number of fragments per particle
  double f_break;
  double volumefraction;
  int n_fragments;

  // stats for breakage
  int n_break,n_break_this,n_break_this_local;
  double mass_break,mass_break_this,mass_break_this_local;

  // data for particles that to be broken
  double **breakdata;
};

}

#endif
#endif
