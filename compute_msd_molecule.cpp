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

#include "compute_msd_molecule.h"
#include "atom.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

ComputeMSDMolecule::ComputeMSDMolecule(LAMMPS *lmp, int narg, char **arg) :
  Compute(lmp, narg, arg)
{
  if (narg != 3) error->all("Illegal compute msd/molecule command");

  if (atom->molecular == 0)
    error->all("Compute msd/molecule requires molecular atom style");

  array_flag = 1;
  size_array_cols = 4;
  extarray = 0;

  // setup molecule-based data and initial COM positions

  nmolecules = molecules_in_group(idlo,idhi);
  size_array_rows = nmolecules;

  massproc = (double *) memory->smalloc(nmolecules*sizeof(double),
					"msd/molecule:massproc");
  masstotal = (double *) memory->smalloc(nmolecules*sizeof(double),
					 "msd/molecule:masstotal");
  com = memory->create_2d_double_array(nmolecules,3,"msd/molecule:com");
  comall = memory->create_2d_double_array(nmolecules,3,"msd/molecule:comall");
  cominit = memory->create_2d_double_array(nmolecules,3,
					   "msd/molecule:cominit");
  msd = memory->create_2d_double_array(nmolecules,4,"msd/molecule:msd");
  array = msd;

  // compute masstotal for each molecule

  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  int i,imol;
  double massone;

  for (i = 0; i < nmolecules; i++) massproc[i] = 0.0;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      massproc[imol] += massone;
    }

  MPI_Allreduce(massproc,masstotal,nmolecules,MPI_DOUBLE,MPI_SUM,world);

  // compute initial COM for each molecule

  firstflag = 1;
  compute_array();
  for (i = 0; i < nmolecules; i++) {
    cominit[i][0] = comall[i][0];
    cominit[i][1] = comall[i][1];
    cominit[i][2] = comall[i][2];
  }
  firstflag = 0;
}

/* ---------------------------------------------------------------------- */

ComputeMSDMolecule::~ComputeMSDMolecule()
{
  memory->sfree(massproc);
  memory->sfree(masstotal);
  memory->destroy_2d_double_array(com);
  memory->destroy_2d_double_array(comall);
  memory->destroy_2d_double_array(cominit);
  memory->destroy_2d_double_array(msd);
}

/* ---------------------------------------------------------------------- */

void ComputeMSDMolecule::init()
{
  int ntmp = molecules_in_group(idlo,idhi);
  if (ntmp != nmolecules)
    error->all("Molecule count changed in compute msd/molecule");
}

/* ---------------------------------------------------------------------- */

void ComputeMSDMolecule::compute_array()
{
  int i,imol;
  double xbox,ybox,zbox,dx,dy,dz;
  double massone;

  invoked_array = update->ntimestep;

  // compute current COM positions

  for (i = 0; i < nmolecules; i++)
    com[i][0] = com[i][1] = com[i][2] = 0.0;

  double **x = atom->x;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int *type = atom->type;
  int *image = atom->image;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;

  double xprd = domain->xprd;
  double yprd = domain->yprd;
  double zprd = domain->zprd;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      xbox = (image[i] & 1023) - 512;
      ybox = (image[i] >> 10 & 1023) - 512;
      zbox = (image[i] >> 20) - 512;
      if (rmass) massone = rmass[i];
      else massone = mass[type[i]];
      imol = molecule[i];
      if (molmap) imol = molmap[imol-idlo];
      else imol--;
      com[imol][0] += (x[i][0] + xbox*xprd) * massone;
      com[imol][1] += (x[i][1] + ybox*yprd) * massone;
      com[imol][2] += (x[i][2] + zbox*zprd) * massone;
    }

  MPI_Allreduce(&com[0][0],&comall[0][0],3*nmolecules,
		MPI_DOUBLE,MPI_SUM,world);
  for (i = 0; i < nmolecules; i++) {
    comall[i][0] /= masstotal[i];
    comall[i][1] /= masstotal[i];
    comall[i][2] /= masstotal[i];
  }

  // MSD is difference between current and initial COM
  // cominit does not yet exist when called from constructor

  if (firstflag) return;

  for (i = 0; i < nmolecules; i++) {
    dx = comall[i][0] - cominit[i][0];
    dy = comall[i][1] - cominit[i][1];
    dz = comall[i][2] - cominit[i][2];
    msd[i][0] = dx*dx;
    msd[i][1] = dy*dy;
    msd[i][2] = dz*dz;
    msd[i][3] = dx*dx + dy*dy + dz*dz;
  }
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeMSDMolecule::memory_usage()
{
  double bytes = 2*nmolecules * sizeof(double);
  if (molmap) bytes += (idhi-idlo+1) * sizeof(int);
  bytes += 2*nmolecules*3 * sizeof(double);
  bytes += nmolecules*4 * sizeof(double);
  return bytes;
}
