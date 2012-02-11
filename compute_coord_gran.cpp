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

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "compute_coord_gran.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "force.h"
#include "pair.h"
#include "comm.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeCoordGran::ComputeCoordGran(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute coord/gran command");

  peratom_flag = 1;
  size_peratom_cols = 0;

  nmax = 0;
  coordination = NULL;

  if(update->ntimestep > 0 && !modify->fix_restart_in_progress()) error->all("Need to define compute coord/gran before first run");
}

/* ---------------------------------------------------------------------- */

ComputeCoordGran::~ComputeCoordGran()
{
  memory->sfree(coordination);
}

/* ---------------------------------------------------------------------- */

void ComputeCoordGran::init()
{
  if (force->pair == NULL)
    error->all("Compute coord/gran requires a pair style be defined");

  if(!atom->radius_flag) error->all("Compute coord/gran requires particles to store radius");

  // need an occasional full neighbor list

  int irequest = neighbor->request((void *) this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->compute = 1;
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->full = 1; 
  neighbor->requests[irequest]->occasional = 1;

  int count = 0;
  for (int i = 0; i < modify->ncompute; i++)
    if (strcmp(modify->compute[i]->style,"coord/gran") == 0) count++;
  if (count > 1 && comm->me == 0)
    error->warning("More than one compute coord/gran");
}

/* ---------------------------------------------------------------------- */

void ComputeCoordGran::init_list(int id, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void ComputeCoordGran::compute_peratom()
{
  int i,j,ii,jj,inum,jnum,n;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq,radsumsq;
  int *ilist,*jlist,*numneigh,**firstneigh;

  invoked_peratom = update->ntimestep;

  // grow coordination array if necessary

  if (atom->nlocal > nmax) {
    memory->sfree(coordination);
    nmax = atom->nmax;
    coordination = (double *)
      memory->smalloc(nmax*sizeof(double),"coord/gran:coordination");
    vector_atom = coordination;
  }

  // invoke neighbor list (will copy or build if necessary)

  neighbor->build_one(list->index);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // compute coordination number for each atom in group

  double **x = atom->x;
  double *radius = atom->radius;
  int *mask = atom->mask;
  int nall = atom->nlocal + atom->nghost;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    if (mask[i] & groupbit) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i];
      jnum = numneigh[i];

      n = 0;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        if (j >= nall) j %= nall;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;
        radsumsq = (radius[i]+radius[j])*(radius[i]+radius[j]);
        if (rsq < radsumsq) n++;
      }

      coordination[i] = n;
    } else coordination[i] = 0.0;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double ComputeCoordGran::memory_usage()
{
  double bytes = nmax * sizeof(double);
  return bytes;
}
