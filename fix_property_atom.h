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

FixStyle(property/atom,FixPropertyAtom)

#else

#ifndef LMP_FIX_PROPERTYPERATOM_H
#define LMP_FIX_PROPERTYPERATOM_H
#include "fix.h"
#include "input.h"

namespace LAMMPS_NS {

enum
{
   FIXPROPERTY_ATOM_SCALAR = 0,
   FIXPROPERTY_ATOM_VECTOR = 1
};

class FixPropertyAtom : public Fix {
 public:
  FixPropertyAtom(class LAMMPS *, int, char **);
  ~FixPropertyAtom();
  int setmask();

  void do_forward_comm();
  void do_reverse_comm();

  Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,bool errflag);

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
  int pack_comm(int, int *, double *, int, int *);
  void unpack_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

  char *variablename;   // name of the variable (used for identification by other fixes)
  int data_style;            // 0 if a scalar is registered, 1 if vector
  int commGhost;        // 1 if communicated to ghost particles (via pack_comm/unpack_comm), 0 if not
  int commGhostRev;     // 1 if rev communicated from ghost particles (via pack_comm_rev/unpack_comm_rev), 0 if not
  int nvalues;
  double *defaultvalues; // default values at particle creation

 private:

}; //end class

}
#endif
#endif
