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

#else

#ifndef LMP_FIX_INSERT_H
#define LMP_FIX_INSERT_H

#include "fix.h"

namespace LAMMPS_NS {

class FixInsert : public Fix {
 public:
  FixInsert(class LAMMPS *, int, char **);
  ~FixInsert();
  virtual int setmask();
  virtual void init();
  virtual void post_create();
  double extend_cut_ghost();
  void pre_exchange();
  virtual void end_of_step() {}

  void write_restart(FILE *);
  virtual void restart(char *);

  virtual double max_rad(int);  
  virtual double max_r_bound();  
  int min_type();
  int max_type();

 protected:

  int iarg;

  /*---TARGET QUANTITIES---how much will be inserted*/

  // insertion target - can bei either number or mass
  int ninsert;
  double massinsert;

  // insertion rate - can bei either number rate or mass rate
  // in particle per time and mass per time units - both can be 0
  double nflowrate;
  double massflowrate;

  //how much we have inserted to far
  int ninserted;
  double massinserted;

  // first, most recent, final insertion step
  int first_ins_step,most_recent_ins_step, final_ins_step;

  /*---INSERTION QUANTITIES---what, where and how exactly will we insert*/

  //particle distribution
  class FixParticledistributionDiscrete *fix_distribution;

  // insert 'ninsert_per' particles every 'insert_every' steps
  int insert_every;
  double ninsert_per;

  // max and min type to be inserted
  int type_max,type_min;

  double maxrad;

  // flag if overlap is checked upon insertion (via all-to-all comm)
  int check_ol_flag;

  // if flag is 1, particles are generated to be in the region as a whole
  // if flag is 0, particles centers are in region
  int all_in_flag;

  int maxattempt;

  // positions generated, and for overlap check
  int nspheres_near;
  double **xnear;

  // velocity and ang vel distribution
  // currently constant - could also be a distribution
  double v_insert[3];
  double omega_insert[3];
  double quat_insert[4];

  /*---FURTHER THINGS THAT WE NEED---*/

  // random generator
  class RanPark *random;
  int seed;

  // all2all comm
  int me,nprocs;
  int *recvcounts,*displs;

  /*---FUNCTION MEMBERS---*/

  virtual void print_stats_start();
  virtual void print_stats_during(int,double);

  virtual void init_defaults();
  virtual void sanity_check();
  virtual void calc_insertion_properties() = 0;

  virtual void pre_insert() {};
  virtual int load_xnear(int);
  virtual int count_nnear();
  virtual int overlap(int) = 0;

  virtual void x_v_omega(int,int&,int&,double&) = 0;

  virtual void finalize_insertion(int){};
};

}

#endif
#endif

