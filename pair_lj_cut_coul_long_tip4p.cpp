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

/* ----------------------------------------------------------------------
   Contributing authors: Amalie Frischknecht and Ahmed Ismail (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_lj_cut_coul_long_tip4p.h"
#include "angle.h"
#include "atom.h"
#include "bond.h"
#include "comm.h"
#include "domain.h"
#include "force.h"
#include "kspace.h"
#include "update.h"
#include "respa.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

#define EWALD_F   1.12837917
#define EWALD_P   0.3275911
#define A1        0.254829592
#define A2       -0.284496736
#define A3        1.421413741
#define A4       -1.453152027
#define A5        1.061405429

/* ---------------------------------------------------------------------- */

PairLJCutCoulLongTIP4P::PairLJCutCoulLongTIP4P(LAMMPS *lmp) : 
  PairLJCutCoulLong(lmp)
{
  single_enable = 0;

  // TIP4P cannot compute virial as F dot r
  // due to find_M() finding bonded H atoms which are not near O atom

  no_virial_compute = 1;
}

/* ---------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::compute(int eflag, int vflag)
{
  int i,j,ii,jj,inum,jnum,itype,jtype,itable;
  int n,vlist[6];
  int iH1,iH2,jH1,jH2;
  double qtmp,xtmp,ytmp,ztmp,delx,dely,delz,evdwl,ecoul;
  double fraction,table;
  double delx1,dely1,delz1,delx2,dely2,delz2,delx3,dely3,delz3;
  double r,r2inv,r6inv,forcecoul,forcelj,cforce,negforce;
  double factor_coul,factor_lj;
  double grij,expm2,prefactor,t,erfc;
  double xiM[3],xjM[3],fO[3],fH[3],v[6];
  double *x1,*x2;
  int *ilist,*jlist,*numneigh,**firstneigh;
  double rsq;

  evdwl = ecoul = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **f = atom->f;
  double **x = atom->x;
  double *q = atom->q;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int nall = nlocal + atom->nghost;
  double *special_coul = force->special_coul;
  double *special_lj = force->special_lj;
  int newton_pair = force->newton_pair;
  double qqrd2e = force->qqrd2e;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  
  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    qtmp = q[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    itype = type[i];
    if (itype == typeO) {
      find_M(i,iH1,iH2,xiM);
      x1 = xiM;
    } else x1 = x[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (j < nall) factor_coul = factor_lj = 1.0;
      else {
	factor_coul = special_coul[j/nall];
	factor_lj = special_lj[j/nall];
	j %= nall;
      }

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      jtype = type[j];
      
      if (rsq < cutsq[itype][jtype]) {

	r2inv = 1.0/rsq;

	if (rsq < cut_ljsq[itype][jtype]) {
	  r6inv = r2inv*r2inv*r2inv;
	  forcelj = r6inv * (lj1[itype][jtype]*r6inv - lj2[itype][jtype]);
	  forcelj *= factor_lj * r2inv;

	  f[i][0] += delx*forcelj;
	  f[i][1] += dely*forcelj;
	  f[i][2] += delz*forcelj;
	  f[j][0] -= delx*forcelj;
	  f[j][1] -= dely*forcelj;
	  f[j][2] -= delz*forcelj;

	  if (eflag) {
	    evdwl = r6inv*(lj3[itype][jtype]*r6inv-lj4[itype][jtype]) -
	      offset[itype][jtype];
	    evdwl *= factor_lj;
	  } else evdwl = 0.0;

	  if (evflag) ev_tally(i,j,nlocal,newton_pair,
			       evdwl,0.0,forcelj,delx,dely,delz);
	}

	// adjust rsq for off-site O charge(s)

	if (itype == typeO || jtype == typeO) { 
	  if (jtype == typeO) {
	    find_M(j,jH1,jH2,xjM);
	    x2 = xjM;
	  } else x2 = x[j];
	  delx = x1[0] - x2[0];
	  dely = x1[1] - x2[1];
	  delz = x1[2] - x2[2];
	  rsq = delx*delx + dely*dely + delz*delz;
	}

	// test current rsq against cutoff and compute Coulombic force

	if (rsq < cut_coulsq) {
	  r2inv = 1 / rsq;
	  if (!ncoultablebits || rsq <= tabinnersq) {
	    r = sqrt(rsq);
	    grij = g_ewald * r;
	    expm2 = exp(-grij*grij);
	    t = 1.0 / (1.0 + EWALD_P*grij);
	    erfc = t * (A1+t*(A2+t*(A3+t*(A4+t*A5)))) * expm2;
	    prefactor = qqrd2e * qtmp*q[j]/r;
	    forcecoul = prefactor * (erfc + EWALD_F*grij*expm2);
	    if (factor_coul < 1.0) {
	      forcecoul -= (1.0-factor_coul)*prefactor; 
	    }
	  } else {
	    union_int_float_t rsq_lookup;
	    rsq_lookup.f = rsq;
	    itable = rsq_lookup.i & ncoulmask;
	    itable >>= ncoulshiftbits;
	    fraction = (rsq_lookup.f - rtable[itable]) * drtable[itable];
	    table = ftable[itable] + fraction*dftable[itable];
	    forcecoul = qtmp*q[j] * table;
	    if (factor_coul < 1.0) {
	      table = ctable[itable] + fraction*dctable[itable];
	      prefactor = qtmp*q[j] * table;
	      forcecoul -= (1.0-factor_coul)*prefactor;
	    }
	  }

	  cforce = forcecoul * r2inv;

	  // if i,j are not O atoms, force is applied directly
	  // if i or j are O atoms, force is on fictitious atom
	  // force partitioning due to Feenstra, J Comp Chem, 20, 786 (1999)
	  // f_f = fictitious force, fO = f_f (1 - 2 alpha), fH = alpha f_f
	  // preserves total force and torque on water molecule
	  // all virial terms are computed as distances to x2
	  // vlist stores 2,4,6 atoms whose forces contribute to virial

	  n = 0;

	  if (itype != typeO) {
	    f[i][0] += delx * cforce;
	    f[i][1] += dely * cforce;
	    f[i][2] += delz * cforce;

	    if (vflag) {
	      v[0] = delx * delx * cforce;
	      v[1] = dely * dely * cforce;
	      v[2] = delz * delz * cforce;
	      v[3] = delx * dely * cforce;
	      v[4] = delx * delz * cforce;
	      v[5] = dely * delz * cforce;

	      vlist[n++] = i;
	    }

	  } else {
	    fO[0] = delx*cforce*(1.0-2.0*alpha);
	    fO[1] = dely*cforce*(1.0-2.0*alpha);
	    fO[2] = delz*cforce*(1.0-2.0*alpha);

	    fH[0] = alpha * (delx*cforce);
	    fH[1] = alpha * (dely*cforce);
	    fH[2] = alpha * (delz*cforce);

	    f[i][0] += fO[0];
	    f[i][1] += fO[1];
	    f[i][2] += fO[2];

	    f[iH1][0] += fH[0];
	    f[iH1][1] += fH[1];
	    f[iH1][2] += fH[2];
	      
	    f[iH2][0] += fH[0];
	    f[iH2][1] += fH[1];
	    f[iH2][2] += fH[2];

	    if (vflag) {
	      delx1 = x[i][0] - x2[0];
	      dely1 = x[i][1] - x2[1];
	      delz1 = x[i][2] - x2[2];
	      domain->minimum_image(delx1,dely1,delz1);

	      delx2 = x[iH1][0] - x2[0];
	      dely2 = x[iH1][1] - x2[1];
	      delz2 = x[iH1][2] - x2[2];
	      domain->minimum_image(delx2,dely2,delz2);

	      delx3 = x[iH2][0] - x2[0];
	      dely3 = x[iH2][1] - x2[1];
	      delz3 = x[iH2][2] - x2[2];
	      domain->minimum_image(delx3,dely3,delz3);

	      v[0] = delx1 * fO[0] + (delx2 + delx3) * fH[0];
	      v[1] = dely1 * fO[1] + (dely2 + dely3) * fH[1];
	      v[2] = delz1 * fO[2] + (delz2 + delz3) * fH[2];
	      v[3] = delx1 * fO[1] + (delx2 + delx3) * fH[1];
	      v[4] = delx1 * fO[2] + (delx2 + delx3) * fH[2];
	      v[5] = dely1 * fO[2] + (dely2 + dely3) * fH[2];

	      vlist[n++] = i;
	      vlist[n++] = iH1;
	      vlist[n++] = iH2;
	    }
	  }

	  if (jtype != typeO) {
	    f[j][0] -= delx * cforce;
	    f[j][1] -= dely * cforce;
	    f[j][2] -= delz * cforce;

	    if (vflag) vlist[n++] = j;

	  } else {
	    negforce = -cforce;

	    fO[0] = delx*negforce*(1.0-2.0*alpha);
	    fO[1] = dely*negforce*(1.0-2.0*alpha);
	    fO[2] = delz*negforce*(1.0-2.0*alpha);

	    fH[0] = alpha * (delx*negforce);
	    fH[1] = alpha * (dely*negforce);
	    fH[2] = alpha * (delz*negforce);

	    f[j][0] += fO[0]; 
	    f[j][1] += fO[1]; 
	    f[j][2] += fO[2]; 
		
	    f[jH1][0] += fH[0];
	    f[jH1][1] += fH[1];
	    f[jH1][2] += fH[2];

	    f[jH2][0] += fH[0];
	    f[jH2][1] += fH[1];
	    f[jH2][2] += fH[2];

	    if (vflag) {
	      delx1 = x[j][0] - x2[0];
	      dely1 = x[j][1] - x2[1];
	      delz1 = x[j][2] - x2[2];
	      domain->minimum_image(delx1,dely1,delz1);

	      delx2 = x[jH1][0] - x2[0];
	      dely2 = x[jH1][1] - x2[1];
	      delz2 = x[jH1][2] - x2[2];
	      domain->minimum_image(delx2,dely2,delz2);

	      delx3 = x[jH2][0] - x2[0];
	      dely3 = x[jH2][1] - x2[1];
	      delz3 = x[jH2][2] - x2[2];
	      domain->minimum_image(delx3,dely3,delz3);

	      v[0] += delx1 * fO[0] + (delx2 + delx3) * fH[0];
	      v[1] += dely1 * fO[1] + (dely2 + dely3) * fH[1];
	      v[2] += delz1 * fO[2] + (delz2 + delz3) * fH[2];
	      v[3] += delx1 * fO[1] + (delx2 + delx3) * fH[1];
	      v[4] += delx1 * fO[2] + (delx2 + delx3) * fH[2];
	      v[5] += dely1 * fO[2] + (dely2 + dely3) * fH[2];

	      vlist[n++] = j;
	      vlist[n++] = jH1;
	      vlist[n++] = jH2;
	    }
	  }

	  if (eflag) {
	    if (!ncoultablebits || rsq <= tabinnersq)
	      ecoul = prefactor*erfc;
	    else {
	      table = etable[itable] + fraction*detable[itable];
	      ecoul = qtmp*q[j] * table;
	    }
	    if (factor_coul < 1.0) ecoul -= (1.0-factor_coul)*prefactor;
	  } else ecoul = 0.0;

	  if (evflag) ev_tally_list(n,vlist,ecoul,v);
	}
      }
    }
  }
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::settings(int narg, char **arg)
{
  if (narg < 6 || narg > 7) error->all("Illegal pair_style command");

  typeO = force->inumeric(arg[0]);
  typeH = force->inumeric(arg[1]);
  typeB = force->inumeric(arg[2]);
  typeA = force->inumeric(arg[3]);
  qdist = force->numeric(arg[4]);

  cut_lj_global = force->numeric(arg[5]);
  if (narg == 6) cut_coul = cut_lj_global;
  else cut_coul = force->numeric(arg[6]);
  
  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
	if (setflag[i][j]) cut_lj[i][j] = cut_lj_global;
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::init_style()
{
  if (atom->tag_enable == 0)
    error->all("Pair style lj/cut/coul/long/tip4p requires atom IDs");
  if (!force->newton_pair) 
    error->all("Pair style lj/cut/coul/long/tip4p requires newton pair on");
  if (!atom->q_flag)
    error->all("Pair style lj/cut/coul/long/tip4p requires atom attribute q");
  if (strcmp(force->kspace_style,"pppm/tip4p") != 0)
    error->all("Pair style is incompatible with KSpace style");
  if (force->bond == NULL)
    error->all("Must use a bond style with TIP4P potential");
  if (force->angle == NULL)
    error->all("Must use an angle style with TIP4P potential");

  PairLJCutCoulLong::init_style();

  // set alpha parameter

  double theta = force->angle->equilibrium_angle(typeA);
  double blen = force->bond->equilibrium_distance(typeB);
  alpha = qdist / (2.0 * cos(0.5*theta) * blen);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::write_restart_settings(FILE *fp)
{
  fwrite(&typeO,sizeof(int),1,fp);
  fwrite(&typeH,sizeof(int),1,fp);
  fwrite(&typeB,sizeof(int),1,fp);
  fwrite(&typeA,sizeof(int),1,fp);
  fwrite(&qdist,sizeof(double),1,fp);

  fwrite(&cut_lj_global,sizeof(double),1,fp);
  fwrite(&cut_coul,sizeof(double),1,fp);
  fwrite(&offset_flag,sizeof(int),1,fp);
  fwrite(&mix_flag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::read_restart_settings(FILE *fp)
{
  if (comm->me == 0) {
    fread(&typeO,sizeof(int),1,fp);
    fread(&typeH,sizeof(int),1,fp);
    fread(&typeB,sizeof(int),1,fp);
    fread(&typeA,sizeof(int),1,fp);
    fread(&qdist,sizeof(double),1,fp);

    fread(&cut_lj_global,sizeof(double),1,fp);
    fread(&cut_coul,sizeof(double),1,fp);
    fread(&offset_flag,sizeof(int),1,fp);
    fread(&mix_flag,sizeof(int),1,fp);
  }

  MPI_Bcast(&typeO,1,MPI_INT,0,world);
  MPI_Bcast(&typeH,1,MPI_INT,0,world);
  MPI_Bcast(&typeB,1,MPI_INT,0,world);
  MPI_Bcast(&typeA,1,MPI_INT,0,world);
  MPI_Bcast(&qdist,1,MPI_DOUBLE,0,world);

  MPI_Bcast(&cut_lj_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&cut_coul,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&offset_flag,1,MPI_INT,0,world);
  MPI_Bcast(&mix_flag,1,MPI_INT,0,world);

}

/* ----------------------------------------------------------------------
  find 2 H atoms bonded to O atom i
  compute position xM of fictitious charge site for O atom
  also return local indices iH1,iH2 of H atoms
------------------------------------------------------------------------- */

void PairLJCutCoulLongTIP4P::find_M(int i, int &iH1, int &iH2, double *xM)
{
  // test that O is correctly bonded to 2 succesive H atoms

  iH1 = atom->map(atom->tag[i] + 1);
  iH2 = atom->map(atom->tag[i] + 2);

  if (iH1 == -1 || iH2 == -1) error->one("TIP4P hydrogen is missing");
  if (atom->type[iH1] != typeH || atom->type[iH2] != typeH)
    error->one("TIP4P hydrogen has incorrect atom type");

  double **x = atom->x; 

  double delx1 = x[iH1][0] - x[i][0];
  double dely1 = x[iH1][1] - x[i][1];
  double delz1 = x[iH1][2] - x[i][2];
  domain->minimum_image(delx1,dely1,delz1);

  double delx2 = x[iH2][0] - x[i][0];
  double dely2 = x[iH2][1] - x[i][1];
  double delz2 = x[iH2][2] - x[i][2];
  domain->minimum_image(delx2,dely2,delz2);

  xM[0] = x[i][0] + alpha * (delx1 + delx2);
  xM[1] = x[i][1] + alpha * (dely1 + dely2);
  xM[2] = x[i][2] + alpha * (delz1 + delz2);
}

/* ---------------------------------------------------------------------- */

void *PairLJCutCoulLongTIP4P::extract(char *str)
{
  if (strcmp(str,"qdist") == 0) return (void *) &qdist;
  else if (strcmp(str,"typeO") == 0) return (void *) &typeO;
  else if (strcmp(str,"typeH") == 0) return (void *) &typeH;
  else if (strcmp(str,"typeA") == 0) return (void *) &typeA;
  else if (strcmp(str,"typeB") == 0) return (void *) &typeB;
  else if (strcmp(str,"cut_coul") == 0) return (void *) &cut_coul;
  return NULL;
}
