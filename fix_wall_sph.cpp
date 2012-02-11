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
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "memory.h"
#include "domain.h"
#include "respa.h"
#include "update.h"
#include "error.h"
#include "sph_kernels.h"
#include "fix_wall_sph.h"

using namespace LAMMPS_NS;

enum{XPLANE,YPLANE,ZPLANE,ZCYLINDER};    // XYZ PLANE need to be 0,1,2

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallSPH::FixWallSPH(LAMMPS *lmp, int narg, char **arg) :
  FixSPH(lmp, narg, arg)
{
  // wallstyle args

  int iarg = 3;
  if (strcmp(arg[iarg],"xplane") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/sph command");
    wallstyle = XPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"yplane") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/sph command");
    wallstyle = YPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zplane") == 0) {
    if (narg < iarg+3) error->all("Illegal fix wall/sph command");
    wallstyle = ZPLANE;
    if (strcmp(arg[iarg+1],"NULL") == 0) lo = -BIG;
    else lo = force->numeric(arg[iarg+1]);
    if (strcmp(arg[iarg+2],"NULL") == 0) hi = BIG;
    else hi = force->numeric(arg[iarg+2]);
    iarg += 3;
  } else if (strcmp(arg[iarg],"zcylinder") == 0) {
    if (narg < iarg+2) error->all("Illegal fix wall/gran command");
    wallstyle = ZCYLINDER;
    lo = hi = 0.0;
    cylradius = force->numeric(arg[iarg+1]);
    iarg += 2;
  }

  // parameters for penetration force
  if (narg < iarg+2) error->all("Illegal fix wall/sph command, not enough arguments for penetration force");
  r0 = force->numeric(arg[iarg]);
  D  = force->numeric(arg[iarg+1]);
  iarg += 2;

  if (wallstyle == XPLANE && domain->xperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == YPLANE && domain->yperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZPLANE && domain->zperiodic)
    error->all("Cannot use wall in periodic dimension");
  if (wallstyle == ZCYLINDER && (domain->xperiodic || domain->yperiodic))
    error->all("Cannot use wall in periodic dimension");

}

/* ---------------------------------------------------------------------- */

FixWallSPH::~FixWallSPH()
{

}

/* ---------------------------------------------------------------------- */

int FixWallSPH::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallSPH::init()
{
  FixSPH::init();

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;
}

/* ---------------------------------------------------------------------- */

void FixWallSPH::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }
}

/* ---------------------------------------------------------------------- */

void FixWallSPH::post_force(int vflag)
{
  double dx,dy,dz,del1,del2,delxy,delr,rsq,r,s,rinv;
  double fwall,gradWmag;

  double wlo = lo;
  double whi = hi;

  double frac,frac2,frac4; // for penetration force

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force on atom if close enough to wall
  //   via wall potential matched to pair potential

  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        if (delxy > 0.) {
          delr = cylradius - delxy;

          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];

        }
      }

      rsq = dx*dx + dy*dy + dz*dz;
      if (rsq == 0.) continue; // center of the cylinder ... no repulsive force!

      r = sqrt(rsq);
      rinv = 1./r;

      // repulsive penetration force
      if (r <= r0) {

        frac = r0*rinv;
        frac2 = frac*frac;
        frac4 = frac2*frac2;

        fwall = D * (frac4 - frac2) * rinv; // second rinv

        f[i][0] += fwall * dx;
        f[i][1] += fwall * dy;
        f[i][2] += fwall * dz;

      }

    }
  } // end loop nlocal
}

/* ---------------------------------------------------------------------- */

void FixWallSPH::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

