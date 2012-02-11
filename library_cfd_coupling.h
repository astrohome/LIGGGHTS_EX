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

/*
   C or Fortran style library interface to LAMMPS
   new LAMMPS-specific functions can be added
*/

#include "mpi.h"

/* ifdefs allow this file to be included in a C program */

#ifdef __cplusplus
extern "C" {
#endif

int liggghts_get_maxtag(void *ptr);
void* locate_coupling_fix(void *ptr);
void data_liggghts_to_of(char *name,char *type,void *ptr,void *&data);
void data_of_to_liggghts(char *name,char *type,void *ptr,void *data);
void update_rm(void *ptr);
void check_datatransfer(void *ptr);

void allocate_external_int(int    **&data, int len2,int len1,int    initvalue,void *ptr);
void allocate_external_double(double **&data, int len2,int len1,double initvalue,void *ptr);

#ifdef __cplusplus
}
#endif

