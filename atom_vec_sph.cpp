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

/* ----------------------------------------------------------------------
Contributing author for SPH:
Andreas Aigner (CD Lab Particulate Flow Modelling, JKU)
andreas.aigner@jku.at
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom_vec_sph.h"
#include "atom.h"
#include "domain.h"
#include "modify.h"
#include "force.h"
#include "fix.h"
#include "fix_adapt.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AtomVecSPH::AtomVecSPH(LAMMPS *lmp, int narg, char **arg) :
  AtomVec(lmp, narg, arg)
{
  molecular = 0;

  comm_x_only = 0; 
  comm_f_only = 1;
  size_forward = 5; 
  size_reverse = 3; 
  size_border = 8;  
  size_velocity = 3; 
  size_data_atom = 6; 
  size_data_vel = 7;
  xcol_data = 5;

  atom->q_flag = atom->density_flag = 1;
  mass_type = 1;

  PI = 4.0*atan(1.0);
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::init()
{
}

/* ----------------------------------------------------------------------
   grow atom arrays
   n = 0 grows arrays by DELTA
   n > 0 allocates arrays to size n
------------------------------------------------------------------------- */

void AtomVecSPH::grow(int n)
{
  if (n == 0) nmax += DELTA;
  else nmax = n;
  atom->nmax = nmax;

  tag = atom->tag = (int *)
    memory->srealloc(atom->tag,nmax*sizeof(int),"atom:tag");
  type = atom->type = (int *)
    memory->srealloc(atom->type,nmax*sizeof(int),"atom:type");
  mask = atom->mask = (int *)
    memory->srealloc(atom->mask,nmax*sizeof(int),"atom:mask");
  image = atom->image = (int *)
    memory->srealloc(atom->image,nmax*sizeof(int),"atom:image");
  x = atom->x = memory->grow_2d_double_array(atom->x,nmax,3,"atom:x");
  v = atom->v = memory->grow_2d_double_array(atom->v,nmax,3,"atom:v");
  f = atom->f = memory->grow_2d_double_array(atom->f,nmax,3,"atom:f");

  q = atom->q = (double *)
    memory->srealloc(atom->q,nmax*sizeof(double),"atom:pressure");
  density = atom->density = (double *)
    memory->srealloc(atom->density,nmax*sizeof(double),"atom:density");

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->grow_arrays(nmax);
}

/* ----------------------------------------------------------------------
   reset local array ptrs
------------------------------------------------------------------------- */

void AtomVecSPH::grow_reset()
{
  tag = atom->tag; type = atom->type;
  mask = atom->mask; image = atom->image;
  x = atom->x; v = atom->v; f = atom->f;
  q = atom->q; density = atom->density;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::copy(int i, int j)
{
  tag[j] = tag[i];
  type[j] = type[i];
  mask[j] = mask[i];
  image[j] = image[i];
  x[j][0] = x[i][0];
  x[j][1] = x[i][1];
  x[j][2] = x[i][2];
  v[j][0] = v[i][0];
  v[j][1] = v[i][1];
  v[j][2] = v[i][2];

  q[j] = q[i];
  density[j] = density[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      modify->fix[atom->extra_grow[iextra]]->copy_arrays(i,j);
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_comm(int n, int *list, double *buf,
			       int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0];
	buf[m++] = x[j][1];
	buf[m++] = x[j][2];
	buf[m++] = q[j];
	buf[m++] = density[j];
      }
    } else {
      if (domain->triclinic == 0) {
	dx = pbc[0]*domain->xprd;
	dy = pbc[1]*domain->yprd;
	dz = pbc[2]*domain->zprd;
      } else {
	dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
	dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
	dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0] + dx;
	buf[m++] = x[j][1] + dy;
	buf[m++] = x[j][2] + dz;
	buf[m++] = q[j];
	buf[m++] = density[j];
      }
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_comm_vel(int n, int *list, double *buf,
				   int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

    m = 0;
    if (pbc_flag == 0) {
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0];
	buf[m++] = x[j][1];
	buf[m++] = x[j][2];
	buf[m++] = q[j];
	buf[m++] = density[j];
	buf[m++] = v[j][0];
	buf[m++] = v[j][1];
	buf[m++] = v[j][2];
      }
    } else {
      if (domain->triclinic == 0) {
	dx = pbc[0]*domain->xprd;
	dy = pbc[1]*domain->yprd;
	dz = pbc[2]*domain->zprd;
      } else {
	dx = pbc[0]*domain->xprd + pbc[5]*domain->xy + pbc[4]*domain->xz;
	dy = pbc[1]*domain->yprd + pbc[3]*domain->yz;
	dz = pbc[2]*domain->zprd;
      }
      for (i = 0; i < n; i++) {
	j = list[i];
	buf[m++] = x[j][0] + dx;
	buf[m++] = x[j][1] + dy;
	buf[m++] = x[j][2] + dz;
	buf[m++] = q[j];
	buf[m++] = density[j];
	buf[m++] = v[j][0];
	buf[m++] = v[j][1];
	buf[m++] = v[j][2];
      }
    }

  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_comm_one(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = density[i];
  return 2;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::unpack_comm(int n, int first, double *buf)
{
  int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      q[i] = buf[m++];
      density[i] = buf[m++];
    }
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::unpack_comm_vel(int n, int first, double *buf)
{
  int i,m,last;

    m = 0;
    last = first + n;
    for (i = first; i < last; i++) {
      x[i][0] = buf[m++];
      x[i][1] = buf[m++];
      x[i][2] = buf[m++];
      q[i] = buf[m++];
      density[i] = buf[m++];
      v[i][0] = buf[m++];
      v[i][1] = buf[m++];
      v[i][2] = buf[m++];
    }

}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::unpack_comm_one(int i, double *buf)
{
  q[i] = buf[0];
  density[i] = buf[1];
  return 2;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_reverse(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = f[i][0];
    buf[m++] = f[i][1];
    buf[m++] = f[i][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::unpack_reverse(int n, int *list, double *buf)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    f[j][0] += buf[m++];
    f[j][1] += buf[m++];
    f[j][2] += buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_border(int n, int *list, double *buf,
				 int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = density[j];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = density[j];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_border_vel(int n, int *list, double *buf,
				     int pbc_flag, int *pbc)
{
  int i,j,m;
  double dx,dy,dz;

  m = 0;
  if (pbc_flag == 0) {
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0];
      buf[m++] = x[j][1];
      buf[m++] = x[j][2];
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = density[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  } else {
    if (domain->triclinic == 0) {
      dx = pbc[0]*domain->xprd;
      dy = pbc[1]*domain->yprd;
      dz = pbc[2]*domain->zprd;
    } else {
      dx = pbc[0];
      dy = pbc[1];
      dz = pbc[2];
    }
    for (i = 0; i < n; i++) {
      j = list[i];
      buf[m++] = x[j][0] + dx;
      buf[m++] = x[j][1] + dy;
      buf[m++] = x[j][2] + dz;
      buf[m++] = tag[j];
      buf[m++] = type[j];
      buf[m++] = mask[j];
      buf[m++] = q[j];
      buf[m++] = density[j];
      buf[m++] = v[j][0];
      buf[m++] = v[j][1];
      buf[m++] = v[j][2];
    }
  }
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::pack_border_one(int i, double *buf)
{
  buf[0] = q[i];
  buf[1] = density[i];
  return 2;
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::unpack_border(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    q[i] = buf[m++];
    density[i] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void AtomVecSPH::unpack_border_vel(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    if (i == nmax) grow(0);
    x[i][0] = buf[m++];
    x[i][1] = buf[m++];
    x[i][2] = buf[m++];
    tag[i] = static_cast<int> (buf[m++]);
    type[i] = static_cast<int> (buf[m++]);
    mask[i] = static_cast<int> (buf[m++]);
    q[i] = buf[m++];
    density[i] = buf[m++];
    v[i][0] = buf[m++];
    v[i][1] = buf[m++];
    v[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::unpack_border_one(int i, double *buf)
{
  q[i] = buf[0];
  density[i] = buf[1];
  return 2;
}

/* ----------------------------------------------------------------------
   pack data for atom I for sending to another proc
   xyz must be 1st 3 values, so comm::exchange() can test on them
------------------------------------------------------------------------- */

int AtomVecSPH::pack_exchange(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];

  buf[m++] = q[i];
  buf[m++] = density[i];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->pack_exchange(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ---------------------------------------------------------------------- */

int AtomVecSPH::unpack_exchange(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);

  q[nlocal] = buf[m++];
  density[nlocal] = buf[m++];

  if (atom->nextra_grow)
    for (int iextra = 0; iextra < atom->nextra_grow; iextra++)
      m += modify->fix[atom->extra_grow[iextra]]->
	unpack_exchange(nlocal,&buf[m]);

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   size of restart data for all atoms owned by this proc
   include extra data stored by fixes
------------------------------------------------------------------------- */

int AtomVecSPH::size_restart()
{
  int i;

  int nlocal = atom->nlocal;
  int n = 11 * nlocal;

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      for (i = 0; i < nlocal; i++)
	n += modify->fix[atom->extra_restart[iextra]]->size_restart(i);

  return n;
}

/* ----------------------------------------------------------------------
   pack atom I's data for restart file including extra quantities
   xyz must be 1st 3 values, so that read_restart can test on them
   molecular types may be negative, but write as positive
------------------------------------------------------------------------- */

int AtomVecSPH::pack_restart(int i, double *buf)
{
  int m = 1;
  buf[m++] = x[i][0];
  buf[m++] = x[i][1];
  buf[m++] = x[i][2];
  buf[m++] = tag[i];
  buf[m++] = type[i];
  buf[m++] = mask[i];
  buf[m++] = image[i];
  buf[m++] = v[i][0];
  buf[m++] = v[i][1];
  buf[m++] = v[i][2];

  if (atom->nextra_restart)
    for (int iextra = 0; iextra < atom->nextra_restart; iextra++)
      m += modify->fix[atom->extra_restart[iextra]]->pack_restart(i,&buf[m]);

  buf[0] = m;
  return m;
}

/* ----------------------------------------------------------------------
   unpack data for one atom from restart file including extra quantities
------------------------------------------------------------------------- */

int AtomVecSPH::unpack_restart(double *buf)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) {
    grow(0);
    if (atom->nextra_store)
      atom->extra = memory->grow_2d_double_array(atom->extra,nmax,
						 atom->nextra_store,
						 "atom:extra");
  }

  int m = 1;
  x[nlocal][0] = buf[m++];
  x[nlocal][1] = buf[m++];
  x[nlocal][2] = buf[m++];
  tag[nlocal] = static_cast<int> (buf[m++]);
  type[nlocal] = static_cast<int> (buf[m++]);
  mask[nlocal] = static_cast<int> (buf[m++]);
  image[nlocal] = static_cast<int> (buf[m++]);
  v[nlocal][0] = buf[m++];
  v[nlocal][1] = buf[m++];
  v[nlocal][2] = buf[m++];

  double **extra = atom->extra;
  if (atom->nextra_store) {
    int size = static_cast<int> (buf[0]) - m;
    for (int i = 0; i < size; i++) extra[nlocal][i] = buf[m++];
  }

  atom->nlocal++;
  return m;
}

/* ----------------------------------------------------------------------
   create one atom of itype at coord
   set other values to defaults
------------------------------------------------------------------------- */

void AtomVecSPH::create_atom(int itype, double *coord)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = 0;
  type[nlocal] = itype;
  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];
  mask[nlocal] = 1;
  image[nlocal] = (512 << 20) | (512 << 10) | 512;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  density[nlocal] = 1.0;
  q[nlocal] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack one line from Atoms section of data file
   initialize other atom quantities
------------------------------------------------------------------------- */

void AtomVecSPH::data_atom(double *coord, int imagetmp, char **values)
{
  int nlocal = atom->nlocal;
  if (nlocal == nmax) grow(0);

  tag[nlocal] = force->inumeric(values[0]);
  if (tag[nlocal] <= 0)
    error->one("Invalid atom ID in Atoms section of data file");

  type[nlocal] = force->inumeric(values[1]);
  if (type[nlocal] <= 0 || type[nlocal] > atom->ntypes)
    error->one("Invalid atom type in Atoms section of data file");

  x[nlocal][0] = coord[0];
  x[nlocal][1] = coord[1];
  x[nlocal][2] = coord[2];

  image[nlocal] = imagetmp;

  mask[nlocal] = 1;
  v[nlocal][0] = 0.0;
  v[nlocal][1] = 0.0;
  v[nlocal][2] = 0.0;

  atom->nlocal++;
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Atoms section of data file
   initialize other atom quantities for this sub-style
------------------------------------------------------------------------- */

int AtomVecSPH::data_atom_hybrid(int nlocal, char **values)
{

  return 0;
}

/* ----------------------------------------------------------------------
   unpack one line from Velocities section of data file
------------------------------------------------------------------------- */

void AtomVecSPH::data_vel(int m, char **values)
{
  v[m][0] = force->numeric(values[0]);
  v[m][1] = force->numeric(values[1]);
  v[m][2] = force->numeric(values[2]);
}

/* ----------------------------------------------------------------------
   unpack hybrid quantities from one line in Velocities section of data file
------------------------------------------------------------------------- */

int AtomVecSPH::data_vel_hybrid(int m, char **values)
{
  return 0;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double AtomVecSPH::memory_usage()
{
  double bytes = 0.0;

  if (atom->memcheck("tag")) bytes += nmax * sizeof(int);
  if (atom->memcheck("type")) bytes += nmax * sizeof(int);
  if (atom->memcheck("mask")) bytes += nmax * sizeof(int);
  if (atom->memcheck("image")) bytes += nmax * sizeof(int);
  if (atom->memcheck("x")) bytes += nmax*3 * sizeof(double);
  if (atom->memcheck("v")) bytes += nmax*3 * sizeof(double);
  if (atom->memcheck("f")) bytes += nmax*3 * sizeof(double);

  if (atom->memcheck("q")) bytes += nmax * sizeof(double);
  if (atom->memcheck("density")) bytes += nmax * sizeof(double);

  return bytes;
}
