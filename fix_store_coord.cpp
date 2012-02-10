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

#include "stdlib.h"
#include "string.h"
#include "fix_store_coord.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixStoreCoord::FixStoreCoord(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 3) error->all("Illegal fix store/coord command");

  restart_peratom = 1;
  peratom_flag = 1;
  size_peratom_cols = 3;
  peratom_freq = 1;

  // optional args

  int comflag = 0;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"com") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix store/coord command");
      if (strcmp(arg[iarg+1],"no") == 0) comflag = 0;
      else if (strcmp(arg[iarg+1],"yes") == 0) comflag = 1;
      else error->all("Illegal fix store/coord command");
      iarg += 2;
    } else error->all("Illegal fix store/coord command");
  }

  // perform initial allocation of atom-based array
  // register with Atom class

  xoriginal = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // cm = original center of mass

  double cm[3];
  if (comflag) {
    double masstotal = group->mass(igroup);
    group->xcm(igroup,masstotal,cm);
  }

  // xoriginal = initial unwrapped positions of atoms
  // relative to center of mass if comflag is set

  double **x = atom->x;
  int *mask = atom->mask;
  int *image = atom->image;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      domain->unmap(x[i],image[i],xoriginal[i]);
      if (comflag) {
	xoriginal[i][0] -= cm[0];
	xoriginal[i][1] -= cm[1];
	xoriginal[i][2] -= cm[2];
      }
    } else xoriginal[i][0] = xoriginal[i][1] = xoriginal[i][2] = 0.0;
  }
}

/* ---------------------------------------------------------------------- */

FixStoreCoord::~FixStoreCoord()
{
  // unregister callbacks to this fix from Atom class
 
  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  memory->destroy_2d_double_array(xoriginal);
}

/* ---------------------------------------------------------------------- */

int FixStoreCoord::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixStoreCoord::memory_usage()
{
  double bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate atom-based array
------------------------------------------------------------------------- */

void FixStoreCoord::grow_arrays(int nmax)
{
  xoriginal =
    memory->grow_2d_double_array(xoriginal,nmax,3,"fix_msd:xoriginal");
  array_atom = xoriginal;
}

/* ----------------------------------------------------------------------
   copy values within local atom-based array
------------------------------------------------------------------------- */

void FixStoreCoord::copy_arrays(int i, int j)
{
  xoriginal[j][0] = xoriginal[i][0];
  xoriginal[j][1] = xoriginal[i][1];
  xoriginal[j][2] = xoriginal[i][2];
}

/* ----------------------------------------------------------------------
   pack values in local atom-based array for exchange with another proc
------------------------------------------------------------------------- */

int FixStoreCoord::pack_exchange(int i, double *buf)
{
  buf[0] = xoriginal[i][0];
  buf[1] = xoriginal[i][1];
  buf[2] = xoriginal[i][2];
  return 3;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based array from exchange with another proc
------------------------------------------------------------------------- */

int FixStoreCoord::unpack_exchange(int nlocal, double *buf)
{
  xoriginal[nlocal][0] = buf[0];
  xoriginal[nlocal][1] = buf[1];
  xoriginal[nlocal][2] = buf[2];
  return 3;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixStoreCoord::pack_restart(int i, double *buf)
{
  buf[0] = 4;
  buf[1] = xoriginal[i][0];
  buf[2] = xoriginal[i][1];
  buf[3] = xoriginal[i][2];
  return 4;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixStoreCoord::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  xoriginal[nlocal][0] = extra[nlocal][m++];
  xoriginal[nlocal][1] = extra[nlocal][m++];
  xoriginal[nlocal][2] = extra[nlocal][m++];
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixStoreCoord::maxsize_restart()
{
  return 4;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixStoreCoord::size_restart(int nlocal)
{
  return 4;
}
