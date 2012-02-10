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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "memory.h"
#include "comm.h"
#include "modify.h"
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling.h"
#include "cfd_datacoupling_mpi.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

CfdDatacouplingMPI::CfdDatacouplingMPI(LAMMPS *lmp,int jarg, int narg, char **arg,FixCfdCoupling* fc) :
  CfdDatacoupling(lmp, jarg, narg, arg,fc)
{
  liggghts_is_active = false;

  if(!atom->tag_enable) error->all("CFD-DEM coupling via MPI requires particles to have tags");

  this->fc = fc;

  nmax_all = 0;
  invoked_last = -1;

  allreduce_short = NULL;
  allreduce_long = NULL;

  if(comm->me == 0) error->info("nevery as specified in LIGGGHTS is overriden by calling external program");

}

CfdDatacouplingMPI::~CfdDatacouplingMPI()
{
    memory->destroy_2d_double_array(allreduce_long);
    memory->destroy_2d_double_array(allreduce_short);
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::pull(char *name,char *type,void *&from)
{
    
    check_grow(type,from);

    int nlocal = atom->nlocal;
    int natomsmax = atom->tag_max();
    int len1 = -1, len2 = -1;

    void * to = find_pull_property(name,type,len1,len2);

    if (!to && nlocal)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s to write data from calling program to.\n",name);
        lmp->error->one("This is fatal");
    }

    if(strcmp(type,"scalar") == 0)
    {
      double **fromdouble = (double**)from;

      if(natomsmax) MPI_Allreduce(&(fromdouble[0][0]),&(allreduce_short[0][0]),natomsmax,MPI_DOUBLE,MPI_SUM,world);

      int m;
      double *todouble = (double*)to;
      for (int i = 0; i < natomsmax; i++) {
        if ((m = atom->map(i+1)) >= 0) {
          todouble[m] = allreduce_short[i][0];
        }
      }
    }

    if(strcmp(type,"vector") == 0)
    {
        double **fromdouble = (double**)from;

        if(natomsmax) MPI_Allreduce(&(fromdouble[0][0]),&(allreduce_long[0][0]),3*natomsmax,MPI_DOUBLE,MPI_SUM,world);

        int m;
        double **todouble = (double**)to;
        for (int i = 0; i < natomsmax; i++) {
          if ((m = atom->map(i+1)) >= 0) {
            todouble[m][0] = allreduce_long[i][0];
            todouble[m][1] = allreduce_long[i][1];
            todouble[m][2] = allreduce_long[i][2];
          }
        }
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::push(char *name,char *type,void *&to)
{

    check_grow(type,to);

    int *tag = atom->tag;
    int nlocal = atom->nlocal;
    int natomsmax = atom->tag_max();

    int len1 = -1, len2 = -1;

    void * from = find_push_property(name,type,len1,len2);

    if (!from && nlocal)
    {
        if(screen) fprintf(screen,"LIGGGHTS could not find property %s requested by calling program.\n",name);
        lmp->error->one("This is fatal");
    }

    int id;
    if(strcmp(type,"scalar") == 0)
    {
        double *fromdouble = (double*)from;
        double **todouble = (double**)to;
        for (int i = 0; i < nlocal; i++) {
            id = tag[i];
            allreduce_short[id-1][0] = fromdouble[i];
        }
        if(natomsmax) MPI_Allreduce(&(allreduce_short[0][0]),&(todouble[0][0]),natomsmax,MPI_DOUBLE,MPI_SUM,world);
    }
    else if(strcmp(type,"vector") == 0)
    {
        double **fromdouble = (double**)from;
        double **todouble = (double**)to;
        for (int i = 0; i < nlocal; i++) {
            id = tag[i];
            allreduce_long[id-1][0] = fromdouble[i][0];
            allreduce_long[id-1][1] = fromdouble[i][1];
            allreduce_long[id-1][2] = fromdouble[i][2];
        }
        if(natomsmax) MPI_Allreduce(&(allreduce_long[0][0]),&(todouble[0][0]),3*natomsmax,MPI_DOUBLE,MPI_SUM,world);
    }
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::check_grow(char *type,void *&ptr)
{
    int tag_max = atom->tag_max();

    if(update->ntimestep > invoked_last)
    {
        
        growflag = false;
        while(nmax_all <= tag_max)
        {
            growflag = true;
            nmax_all += DELTA;
        }

        if(growflag)
        {
            grow_arrays_allred(nmax_all);
            
        }

        invoked_last = update->ntimestep;

        reset_arrays_allred(nmax_all);
    }

    else if(growflag)
    {
        
    }
}

/* ---------------------------------------------------------------------- */

void CfdDatacouplingMPI::allocate_external(int **&data, int len2,int len1,int initvalue)
{
  if(len1 == -1) len1 = atom->tag_max();
  data = (int**)(lmp->memory->grow_2d_int_array(data, len1,len2, "CfdDatacouplingMPI:data"));
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

void CfdDatacouplingMPI::allocate_external(double **&data, int len2,int len1,double initvalue)
{
  if(len1 == -1) len1 = atom->tag_max();
  data = (double**)(lmp->memory->grow_2d_double_array(data, len1,len2, "CfdDatacouplingMPI:data"));
  for (int i = 0; i < len1; i++)
    for (int j = 0; j < len2; j++)
        data[i][j] = initvalue;
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::grow_arrays_allred(int nmax)
{
  allreduce_short = (double**)(lmp->memory->grow_2d_double_array(allreduce_short, nmax,1, "fix_cfd_coupling:allreduce_short"));
  allreduce_long = (double**)(lmp->memory->grow_2d_double_array(allreduce_long, nmax,3, "fix_cfd_coupling:allreduce_long"));
  reset_arrays_allred(nmax);
}

/* ---------------------------------------------------------------------- */

inline void CfdDatacouplingMPI::reset_arrays_allred(int nmax)
{
  
  for (int ii = 0; ii < nmax; ii++)
      allreduce_short[ii][0] = 0.;
  for (int ii = 0; ii < nmax; ii++)
    for (int jj = 0; jj < 3; jj++)
      allreduce_long[ii][jj]=0.;
}

