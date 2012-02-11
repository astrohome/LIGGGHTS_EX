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

FixStyle(property/global,FixPropertyGlobal)
FixStyle(property/atomtype,FixPropertyGlobal)
FixStyle(property/atomtypepair,FixPropertyGlobal)

#else

#ifndef LMP_FIX_PROPERTYGLOBAL_H
#define LMP_FIX_PROPERTYGLOBAL_H
#include "fix.h"
#include "input.h"

namespace LAMMPS_NS {

enum
{
        FIXPROPERTY_GLOBAL_SCALAR = 0,
        FIXPROPERTY_GLOBAL_VECTOR = 1,
        FIXPROPERTY_GLOBAL_MATRIX = 2
};

class FixPropertyGlobal : public Fix {
 friend class Modify;
 friend class MechParamGran;
 friend class CfdDatacouplingFile;

 public:
  FixPropertyGlobal(class LAMMPS *, int, char **);
  ~FixPropertyGlobal();
  int setmask();
  void init();

  Fix* check_fix(const char *varname,const char *svmstyle,int len1,int len2,bool errflag);

  double memory_usage();
  double compute_scalar();
  double compute_vector(int);
  double compute_vector_modified(int);
  double compute_array(int,int);
  double compute_array_modified(int,int);
  void vector_modify(int,double);
  void array_modify(int,int,double);
  void new_array(int l1,int l2);

  //bool checkCorrectness(int,char*,int,int);

  const double* get_values() {return values;}
  const double* get_values_modified() {return values_recomputed;}
  double const* const* get_array() {return array;}
  double const* const* get_array_modified() {return array_recomputed;}

  void grow(int,int);

  char *variablename;   // name of the variable (used for identification by other fixes)
  int data_style;         // 0 if a scalar is registered, 1 if vector, 2 if 2d array (matrix)
  int nvalues;
  double *values;            // original values to be stored in this fix
  double *values_recomputed; // values to be stored in this fix, recomputed by eg another fix
  double **array;
  double **array_recomputed;

}; //end class

}
#endif
#endif

