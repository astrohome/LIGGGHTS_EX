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

#ifndef LMP_FIX_CFD_COUPLING_H
#define LMP_FIX_CFD_COUPLING_H

#include "fix.h"

namespace LAMMPS_NS {

class FixCfdCoupling : public Fix {
 friend class CfdRegionmodel;
 friend class CfdDatacoupling;
 public:
  FixCfdCoupling(class LAMMPS *, int, char **);
  ~FixCfdCoupling();

  virtual int setmask();
  void init();
  virtual void init_submodel(){}
  virtual void setup(int);
  virtual void min_setup(int);
  virtual void initial_integrate(int) {}
  virtual void post_integrate() {}
  virtual void pre_exchange() {}
  virtual void pre_neighbor() {}
  virtual void pre_force(int) {}
  virtual void post_force(int) {}
  virtual void final_integrate() {}
  void end_of_step();
  void post_force_respa(int, int, int);
  void min_post_force(int);

  void pull(char *name,char *type,void *&ptr);
  void push(char *name,char *type,void *&ptr);

  void add_push_property(char *name,char *type);
  void add_push_prop(char *name,char *type);
  void add_pull_property(char *name,char *type);
  void add_pull_prop(char *name,char *type);

  virtual void special_settings() =0;
  virtual void post_create();
  virtual void pre_delete(){}

  int coupleThis() {return master->couple_this;}
  bool is_master() {return master==this;}

  void check_datatransfer();
  void allocate_external(int    **&data, int len2,int len1,int    initvalue);
  void allocate_external(double **&data, int len2,int len1,double initvalue);

 protected:
  int couple_this;
  class FixCfdCoupling *master;
  int master_flag;
  int iarg;

 protected:

  void* find_property(int flag,char *name,char *type,int &len1,int &len2);
  void* find_pull_property(char *name,char *type,int &len1,int &len2);
  void* find_push_property(char *name,char *type,int &len1,int &len2);

  int couple_nevery,ts_create;

  int nvalues_max;

  int npull;
  
  char **pullnames;
  char **pulltypes;
  
  int *pullinvoked;
  
  int npush;
  char **pushnames;
  char **pushtypes;
  
  int *pushinvoked;
  
  void grow_();

  class CfdDatacoupling *dc;

  class CfdRegionmodel *rm;

  int nlevels_respa;
};

}

#endif
#endif
