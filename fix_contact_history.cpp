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
#include "stdio.h"
#include "stdlib.h"
#include "fix_contact_history.h"
#include "atom.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "pair_gran.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define DELTA_MAXTOUCH 15

/* ---------------------------------------------------------------------- */

FixContactHistory::FixContactHistory(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  //parse args
  if(narg < 4) error->all("Illegal fix contacthistory command - not enough parameters");

  //read dnum
  dnum = atoi(arg[3]);
  
  if(dnum < 1) error->all("dnum must be >=1 in fix contacthistory");
  if(dnum > 10) error->warning("dnum >10 in fix contacthistory - are you really sure you intend this?");

  //read newtonflag
  if(narg-4 < 2*dnum) error->all("Illegal fix contacthistory command - not enough parameters (need to specify an id and a newtonflag for each dnum)");

  int iarg = 4;
  newtonflag = new int[dnum];
  //history_id = new char*[dnum];
  history_id = (char**) memory->smalloc((dnum)*sizeof(char*),"FixContactHistory:history_id");
  for(int i = 0 ; i < dnum; i++)
  {
    
    history_id[i] = new char[strlen(arg[iarg])+1];
    strcpy(history_id[i],arg[iarg++]);
    newtonflag[i] = atoi(arg[iarg++]);
    if(newtonflag[i] != 0 && newtonflag[i] != 1) error->all("Illegal fix history command - newtonflag must be either 0 or 1");
  }

  restart_peratom = 1;
  restart_global = 1; 
  create_attribute = 1;

  // set time_depend so that history will be preserved correctly
  // across multiple runs via laststep setting in granular pair styles
  time_depend = 1;

  maxtouch = DELTA_MAXTOUCH;

  // perform initial allocation of atom-based arrays
  // register with atom class

  npartner = NULL;
  partner = NULL;
  contacthistory = NULL;
  grow_arrays(atom->nmax);
  atom->add_callback(0);
  atom->add_callback(1);

  // initialize npartner to 0 so neighbor list creation is OK the 1st time

  int nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; i++) npartner[i] = 0;

  laststep = -1;
}

/* ---------------------------------------------------------------------- */

FixContactHistory::~FixContactHistory()
{
  // unregister this fix so atom class doesn't invoke it any more

  atom->delete_callback(id,0);
  atom->delete_callback(id,1);

  // delete locally stored arrays

  memory->sfree(npartner);
  memory->destroy_2d_int_array(partner);
  memory->destroy_3d_double_array(contacthistory);

  delete[]newtonflag;

  for(int i = 0; i < dnum; i++) delete [](history_id[i]);
  memory->sfree(history_id);
}

/* ---------------------------------------------------------------------- */

int FixContactHistory::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::init()
{
    if (atom->tag_enable == 0)
      error->all("Pair style granular with history requires atoms have IDs");

    if(!force->pair_match("gran", 0)) error->all("Please use a granular pair style for fix contacthistory");
    pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
}

/* ---------------------------------------------------------------------- */

void FixContactHistory::setup_pre_exchange()
{
  if(laststep > 0) pre_exchange();
}

/* ----------------------------------------------------------------------
   copy shear partner info from neighbor lists to atom arrays
   so can be exchanged with atoms
------------------------------------------------------------------------- */

void FixContactHistory::pre_exchange()
{
  int i,j,ii,jj,m,inum,jnum,d;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *hist,*allhist,**firsthist;

  // zero npartners for all current atoms

  int nlocal = atom->nlocal;
  for (i = 0; i < nlocal; i++) npartner[i] = 0;

  // copy contact info from neighbor list atoms to atom arrays

  int *tag = atom->tag;

  NeighList *list = pair_gran->list;
  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = list->listgranhistory->firstneigh;
  firsthist = list->listgranhistory->firstdouble;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    jlist = firstneigh[i];
    allhist = firsthist[i];
    jnum = numneigh[i];
    touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      if (touch[jj]) {
        hist = &allhist[dnum*jj]; 
        j = jlist[jj];

        if (npartner[i] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 
        m = npartner[i];
        partner[i][m] = tag[j];
        for (d = 0; d < dnum; d++) {
            contacthistory[i][m][d] = hist[d]; 
        }
        npartner[i]++;
        if (j < nlocal) {
            if (npartner[j] >= maxtouch) grow_arrays_maxtouch(atom->nmax); 

            m = npartner[j];
            partner[j][m] = tag[i];
            for (d = 0; d < dnum; d++) {
               if(newtonflag[d]) contacthistory[j][m][d] = -hist[d]; 
               else              contacthistory[j][m][d] =  hist[d]; 
            }
          npartner[j]++;
        }
      }
    }
  }

  int maxtouch_all;
  MPI_Allreduce(&maxtouch,&maxtouch_all,1,MPI_INT,MPI_MAX,world);
  while (maxtouch<maxtouch_all) grow_arrays_maxtouch(atom->nmax);

  laststep = update->ntimestep;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixContactHistory::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = nmax * sizeof(int);
  bytes += nmax*maxtouch * sizeof(int);  
  bytes += nmax*maxtouch * dnum * sizeof(double); 
  return bytes;
}

/* ----------------------------------------------------------------------
   allocate local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistory::grow_arrays(int nmax)
{
  npartner = (int *) memory->srealloc(npartner,nmax*sizeof(int),
				      "contacthistory:npartner");

  partner = memory->grow_2d_int_array(partner,nmax,maxtouch,
				      "contacthistory:partner");
  //fprintf(screen,"nmax=%d, maxtouch=%d\n",nmax,maxtouch);
  contacthistory =
    memory->grow_3d_double_array(contacthistory,nmax,maxtouch,dnum,
				 "contact_history:contacthistory");
}

/* ----------------------------------------------------------------------
   grow local atom-based arrays in case maxtouch is too small 
------------------------------------------------------------------------- */

void FixContactHistory::grow_arrays_maxtouch(int nmax)
{
  if(comm->me==0)
  {
      if(screen) fprintf(screen,"INFO: more than %d touching neighbor atoms found, growing contact history.\n",maxtouch);
      if(logfile) fprintf(logfile,"INFO: more than %d touching neighbor atoms found, growing contact history.\n",maxtouch);
  }

  int **partner_g = memory->create_2d_int_array(nmax,maxtouch+DELTA_MAXTOUCH,
				      "contacthistory:partner_g");
  double ***contacthistory_g =
    memory->create_3d_double_array(nmax,maxtouch+DELTA_MAXTOUCH,dnum,"shear_history:contacthistory_g");

  for (int i = 0; i < nmax; i++)
  {
      for (int j = 0; j < maxtouch; j++)
      {
          partner_g[i][j] = partner[i][j];
          for (int k = 0 ; k < dnum; k++)
            contacthistory_g[i][j][k] = contacthistory_g[i][j][k];
      }
  }
  maxtouch += DELTA_MAXTOUCH;
  int **h1; double ***h2;			 ;
  h1 = partner;
  h2 = contacthistory;
  partner = partner_g;
  contacthistory = contacthistory_g;
  memory->destroy_2d_int_array(h1);
  memory->destroy_3d_double_array(h2);
}

/* ----------------------------------------------------------------------
   copy values within local atom-based arrays
------------------------------------------------------------------------- */

void FixContactHistory::copy_arrays(int i, int j)
{
  npartner[j] = npartner[i];
  for (int m = 0; m < npartner[j]; m++) {
    partner[j][m] = partner[i][m];
    for (int d = 0; d < dnum; d++) {
      contacthistory[j][m][d] = contacthistory[i][m][d];
    }
  }
}

/* ----------------------------------------------------------------------
   initialize one atom's array values, called when atom is created
------------------------------------------------------------------------- */

void FixContactHistory::set_arrays(int i)
{
  npartner[i] = 0;
}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    for (int d = 0; d < dnum; d++) {
      buf[m++] = contacthistory[i][n][d];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values in local atom-based arrays from exchange with another proc
------------------------------------------------------------------------- */

int FixContactHistory::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  npartner[nlocal] = static_cast<int> (buf[m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (buf[m++]);
    for (int d = 0; d < dnum; d++) {
      contacthistory[nlocal][n][d] = buf[m++];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixContactHistory::write_restart(FILE *fp)
{
  int n = 0;
  double list[6];
  list[n++] = static_cast<double>(dnum);

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }

  pre_exchange();
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixContactHistory::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  int pack_dnum = static_cast<int> (list[n++]);
  if(pack_dnum != dnum) error->all("Saved simulation state used different contact history model - can not restart");

}

/* ----------------------------------------------------------------------
   pack values in local atom-based arrays for restart file
------------------------------------------------------------------------- */

int FixContactHistory::pack_restart(int i, double *buf)
{
  int m = 0;
  buf[m++] = (dnum+1)*npartner[i] + 2;
  buf[m++] = npartner[i];
  for (int n = 0; n < npartner[i]; n++) {
    buf[m++] = partner[i][n];
    for (int d = 0; d < dnum; d++) {
      buf[m++] = contacthistory[i][n][d];
    }
  }
  return m;
}

/* ----------------------------------------------------------------------
   unpack values from atom->extra array to restart the fix
------------------------------------------------------------------------- */

void FixContactHistory::unpack_restart(int nlocal, int nth)
{
  double **extra = atom->extra;

  // skip to Nth set of extra values

  int m = 0;
  for (int i = 0; i < nth; i++) m += static_cast<int> (extra[nlocal][m]);
  m++;

  int d;

  npartner[nlocal] = static_cast<int> (extra[nlocal][m++]);
  for (int n = 0; n < npartner[nlocal]; n++) {
    partner[nlocal][n] = static_cast<int> (extra[nlocal][m++]);
    for (d = 0; d < dnum; d++) {
      contacthistory[nlocal][n][d] = extra[nlocal][m++];
    }
  }
}

/* ----------------------------------------------------------------------
   maxsize of any atom's restart data
------------------------------------------------------------------------- */

int FixContactHistory::maxsize_restart()
{
  return (dnum+1)*maxtouch + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixContactHistory::size_restart(int nlocal)
{
  return (dnum+1)*npartner[nlocal] + 2;
}

/* ----------------------------------------------------------------------
   size of atom nlocal's restart data
------------------------------------------------------------------------- */

int FixContactHistory::n_contacts()
{
    int ncontacts = 0;
    int nlocal = atom->nlocal;

    for(int i = 0; i < nlocal; i++)
        ncontacts += npartner[i];

    return ncontacts;
}
