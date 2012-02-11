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

#ifndef LMP_INPUT_MULTISPHERE_H
#define LMP_INPUT_MULTISPHERE_H

#include "stdio.h"
#include "input.h"

namespace LAMMPS_NS {

class InputMultisphere : protected Input {
 public:

  InputMultisphere(class LAMMPS *, int, char **);
  ~InputMultisphere();

  int clmpfile(double **,double*,int);
  void clmpfile(const char *,double **,double*,int);

};
}

#endif

