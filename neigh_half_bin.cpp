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

#include "neighbor.h"
#include "neigh_list.h"
#include "atom.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   binned neighbor list construction with partial Newton's 3rd law
   each owned atom i checks own bin and other bins in stencil
   pair stored once if i,j are both owned and i < j
   pair stored by me if j is ghost (also stored by proc owning j)
------------------------------------------------------------------------- */

void Neighbor::half_bin_no_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  if (includegroup) nlocal = atom->nfirst;
  int molecular = atom->molecular;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in other bins in stencil including self
    // only store pair if i < j
    // stores own/own pairs only once
    // stores own/ghost pairs on both procs

    ibin = coord2bin(x[i]);

    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
        if (j <= i) continue;

        jtype = type[j];
        if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

        delx = xtmp - x[j][0];
        dely = ytmp - x[j][1];
        delz = ztmp - x[j][2];
        rsq = delx*delx + dely*dely + delz*delz;

        if (rsq <= cutneighsq[itype][jtype]) {
          if (molecular) which = find_special(i,j);
          else which = 0;
          if (which == 0) neighptr[n++] = j;
          else if (which > 0) neighptr[n++] = which*nall + j;
        }
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with full Newton's 3rd law
   each owned atom i checks its own bin and other bins in Newton stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int molecular = atom->molecular;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over rest of atoms in i's bin, ghosts are at end of linked list
    // if j is owned atom, store it, since j is beyond i in linked list
    // if j is ghost, only store if j coords are "above and to the right" of i

    for (j = bins[i]; j >= 0; j = bins[j]) {
      if (j >= nlocal) {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] < xtmp) continue;
      }

      jtype = type[j];
      if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;

      if (rsq <= cutneighsq[itype][jtype]) {
	if (molecular) which = find_special(i,j);
	else which = 0;
	if (which == 0) neighptr[n++] = j;
	else if (which > 0) neighptr[n++] = which*nall + j;
      }
    }

    // loop over all atoms in other bins in stencil, store every pair

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	jtype = type[j];
	if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;
	}
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}

/* ----------------------------------------------------------------------
   binned neighbor list construction with Newton's 3rd law for triclinic
   each owned atom i checks its own bin and other bins in triclinic stencil
   every pair stored exactly once by some processor
------------------------------------------------------------------------- */

void Neighbor::half_bin_newton_tri(NeighList *list)
{
  int i,j,k,n,itype,jtype,ibin,which;
  double xtmp,ytmp,ztmp,delx,dely,delz,rsq;
  int *neighptr;

  // bin local & ghost atoms

  bin_atoms();

  // loop over each atom, storing neighbors

  double **x = atom->x;
  int *type = atom->type;
  int *mask = atom->mask;
  int *molecule = atom->molecule;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  int molecular = atom->molecular;
  if (includegroup) nlocal = atom->nfirst;

  int *ilist = list->ilist;
  int *numneigh = list->numneigh;
  int **firstneigh = list->firstneigh;
  int **pages = list->pages;
  int nstencil = list->nstencil;
  int *stencil = list->stencil;

  int inum = 0;
  int npage = 0;
  int npnt = 0;

  for (i = 0; i < nlocal; i++) {

    if (pgsize - npnt < oneatom) {
      npnt = 0;
      npage++;
      if (npage == list->maxpage) pages = list->add_pages();
    }

    neighptr = &pages[npage][npnt];
    n = 0;

    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // loop over all atoms in bins in stencil
    // pairs for atoms j "below" i are excluded
    // below = lower z or (equal z and lower y) or (equal zy and <= x)
    // this excludes self-self interaction

    ibin = coord2bin(x[i]);
    for (k = 0; k < nstencil; k++) {
      for (j = binhead[ibin+stencil[k]]; j >= 0; j = bins[j]) {
	if (x[j][2] < ztmp) continue;
	if (x[j][2] == ztmp && x[j][1] < ytmp) continue;
	if (x[j][2] == ztmp && x[j][1] == ytmp && x[j][0] <= xtmp) continue;

	jtype = type[j];
	if (exclude && exclusion(i,j,itype,jtype,mask,molecule)) continue;

	delx = xtmp - x[j][0];
	dely = ytmp - x[j][1];
	delz = ztmp - x[j][2];
	rsq = delx*delx + dely*dely + delz*delz;

	if (rsq <= cutneighsq[itype][jtype]) {
	  if (molecular) which = find_special(i,j);
	  else which = 0;
	  if (which == 0) neighptr[n++] = j;
	  else if (which > 0) neighptr[n++] = which*nall + j;
	}
      }
    }

    ilist[inum++] = i;
    firstneigh[i] = neighptr;
    numneigh[i] = n;
    npnt += n;
    if (n > oneatom || npnt >= pgsize)
      error->one("Neighbor list overflow, boost neigh_modify one or page");
  }

  list->inum = inum;
}
