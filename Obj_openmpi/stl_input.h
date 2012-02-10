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

#ifndef LMP_INPUT_STL_H
#define LMP_INPUT_STL_H

#include "stdio.h"
#include "input.h"

namespace LAMMPS_NS {

class InputSTL : protected Input {
 public:

  InputSTL(class LAMMPS *, int, char **);
  ~InputSTL();

  void stlfile(class FixMeshGran *);                    // process stl file modified C.K.
  void stlfile(const char *,class FixMeshGran *);       // process stl file modified C.K.
};

}

#endif

