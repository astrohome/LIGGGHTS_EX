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

#ifdef COMPUTE_CLASS

ComputeStyle(coord/gran,ComputeCoordGran)

#else

#ifndef LMP_COMPUTE_COORD_GRAN_H
#define LMP_COMPUTE_COORD_GRAN_H

#include "compute.h"

namespace LAMMPS_NS {

class ComputeCoordGran : public Compute {
 public:
  ComputeCoordGran(class LAMMPS *, int, char **);
  ~ComputeCoordGran();
  void init();
  void init_list(int, class NeighList *);
  void compute_peratom();
  double memory_usage();

 private:
  int nmax;
  class NeighList *list;
  double *coordination;
};

}

#endif
#endif
