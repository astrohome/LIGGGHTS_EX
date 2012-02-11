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

#ifdef COMMAND_CLASS

CommandStyle(read_restart,ReadRestart)

#else

#ifndef LMP_READ_RESTART_H
#define LMP_READ_RESTART_H

#include "stdio.h"
#include "pointers.h"

namespace LAMMPS_NS {

class ReadRestart : protected Pointers {
 public:
  ReadRestart(class LAMMPS *);
  void command(int, char **);

 private:
  int me;
  FILE *fp;
  int nprocs_file;

  void file_search(char *, char *);
  void header();
  void type_arrays();
  void force_fields();

  int read_int();
  double read_double();
  char *read_char();
};

}

#endif
#endif
