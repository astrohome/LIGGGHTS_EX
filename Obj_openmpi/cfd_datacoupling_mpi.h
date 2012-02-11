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

#ifdef CFD_DATACOUPLING_CLASS

   CfdDataCouplingStyle(mpi,CfdDatacouplingMPI)

#else

#ifndef LMP_CFD_DATACOUPLING_MPI_H
#define LMP_CFD_DATACOUPLING_MPI_H

#include "cfd_datacoupling.h"

namespace LAMMPS_NS {

class CfdDatacouplingMPI : public CfdDatacoupling {
 public:
  CfdDatacouplingMPI(class LAMMPS *, int,int, char **,class FixCfdCoupling*);
  ~CfdDatacouplingMPI();

  void pull(char *,char *,void *&);
  void push(char *,char *,void *&);

  void allocate_external(int    **&data, int len2,int len1,int    initvalue);
  void allocate_external(double **&data, int len2,int len1,double initvalue);

 private:
  
  int invoked_last;
  int nmax_all;
  bool growflag;
  double **allreduce_long,**allreduce_short;

  void check_grow(char *type,void *&ptr);
  void grow_arrays_allred(int nmax);
  void reset_arrays_allred(int nmax);
};

}

#endif
#endif
