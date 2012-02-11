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

FixStyle(contacthistory,FixContactHistory)

#else

#ifndef LMP_FIX_CONTACT_HISTORY_H
#define LMP_FIX_CONTACT_HISTORY_H

#include "fix.h"

namespace LAMMPS_NS {

class FixContactHistory : public Fix {
  friend class Neighbor;
  friend class PairGran;
  friend class CMGranModelframework;

 public:
  FixContactHistory(class LAMMPS *, int, char **);
  ~FixContactHistory();
  int setmask();
  void init();
  void setup_pre_exchange();
  void pre_exchange();

  double memory_usage();
  void grow_arrays(int);
  void copy_arrays(int, int);
  void set_arrays(int);
  int pack_exchange(int, double *);
  int unpack_exchange(int, double *);
  int pack_restart(int, double *);
  void unpack_restart(int, int);
  int size_restart(int);
  int maxsize_restart();
  void write_restart(FILE *);
  void restart(char *);

  int n_contacts(); 

 private:
  int *npartner;                // # of touching partners of each atom
  int **partner;                // tags for the partners
  double ***contacthistory;     // history values with the partner
  int maxtouch;                 // max number of partners per atom 
  void grow_arrays_maxtouch(int);

  class PairGran *pair_gran;

  int dnum;
  int *newtonflag;
  char **history_id;
  int laststep;
};

}

#endif
#endif
