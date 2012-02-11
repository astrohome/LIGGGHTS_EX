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
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_pour.h"
#include "fix_contact_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "mech_param_gran.h"
#include "fix_rigid.h"
#include "fix_pour.h"
#include "fix_pour_dev.h"
#include "fix_particledistribution_discrete.h"
#include "fix_pour_legacy.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "compute_pair_gran_local.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGran::PairGran(LAMMPS *lmp) : Pair(lmp)
{
  single_enable = 0;
  no_virial_compute = 1;

  mpg = new MechParamGran(lmp);

  history = 0;
  fix_history = NULL;

  cpl_enable = 1;
  cpl = NULL;

  energytrack_enable = 0;
  fppaCPEn = fppaCDEn = fppaCPEt = fppaCDEVt = fppaCDEFt = fppaCTFW = fppaDEH = NULL;
  CPEn = CDEn = CPEt = CDEVt = CDEFt = CTFW = DEH = NULL;

  laststep = -1;
}

/* ---------------------------------------------------------------------- */

PairGran::~PairGran()
{
  if (fix_history) modify->delete_fix("contacthistory");

  if (allocated) {
    memory->destroy_2d_int_array(setflag);
    memory->destroy_2d_double_array(cutsq);
  }
  delete mpg;

  //unregister energy terms as property/atom
  if (fppaCPEn) modify->delete_fix("CPEn");
  if (fppaCDEn) modify->delete_fix("CDEn");
  if (fppaCPEt) modify->delete_fix("CPEt");
  if (fppaCDEVt) modify->delete_fix("CDEVt");
  if (fppaCDEFt) modify->delete_fix("CDEFt");
  if (fppaCTFW) modify->delete_fix("CTFW");
  if (fppaDEH) modify->delete_fix("DEH");

  delete [] onerad_dynamic;
  delete [] onerad_frozen;
  delete [] maxrad_dynamic;
  delete [] maxrad_frozen;
}

/* ---------------------------------------------------------------------- */

void PairGran::updatePtrs()
{
	if(fppaCPEn) CPEn = fppaCPEn->vector_atom;
	if(fppaCDEn) CDEn = fppaCDEn->vector_atom;
	if(fppaCPEt) CPEt = fppaCPEt->vector_atom;
	if(fppaCDEVt) CDEVt = fppaCDEVt->vector_atom;
	if(fppaCDEFt) CDEFt = fppaCDEFt->vector_atom;
	if(fppaCTFW) CTFW = fppaCTFW->vector_atom;
	if(fppaDEH) DEH = fppaDEH->vector_atom;
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairGran::allocate()
{
  allocated = 1;
  int n = atom->ntypes; //ensured elsewhere that this is high enough

  setflag = memory->create_2d_int_array(n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  cutsq = memory->create_2d_double_array(n+1,n+1,"pair:cutsq");

  onerad_dynamic = new double[n+1];
  onerad_frozen = new double[n+1];
  maxrad_dynamic = new double[n+1];
  maxrad_frozen = new double[n+1];
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairGran::coeff(int narg, char **arg)
{
  if (narg > 2) error->all("Incorrect args for pair coefficients");
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi;
  force->bounds(arg[0],atom->ntypes,ilo,ihi);
  force->bounds(arg[1],atom->ntypes,jlo,jhi);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    for (int j = MAX(jlo,i); j <= jhi; j++) {
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all("Incorrect args for pair coefficients");
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairGran::init_style()
{
  int i, nRigid = 0;

  // error and warning checks

  if(strcmp(update->unit_style,"metal")==0 || strcmp(update->unit_style,"real")==0) error->all("Cannot use a non-consistent unit system with pair gran. Please use si,cgs or lj.");

  if (!atom->radius_flag || !atom->omega_flag || !atom->torque_flag)
    error->all("Pair granular requires atom attributes radius, omega, torque");
  if (comm->ghost_velocity == 0)
    error->all("Pair granular requires ghost atoms store velocity");

  // check if a fix rigid is registered - important for damp
  
  fr = NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
  {
      if(strcmp(modify->fix[ifix]->style,"rigid") == 0)
      {
          fr = static_cast<FixRigid*>(modify->fix[ifix]);
          nRigid++;
      }
  }

  if(nRigid > 1) error->warning("Pair gran does currently not support more than one fix rigid. This may result in under-damping.");

  dt = update->dt;

  // if shear history is stored:
  // check if newton flag is valid
  // if first init, create Fix needed for storing shear history

  if (history && force->newton_pair == 1)
    error->all("Pair granular with shear history requires newton pair off");

  // init contact history
  if(history)
  {
    if (!fix_history)
    {
        
        char **fixarg = new char*[4+2*dnum];
        fixarg[3] = new char[3];
        history_args(&fixarg[4]);

        fixarg[0] = (char *) "contacthistory";
        fixarg[1] = (char *) "all";
        fixarg[2] = (char *) "contacthistory";
        sprintf(fixarg[3],"%d",dnum);
        modify->add_fix(4+2*dnum,fixarg);

        if(modify->n_fixes_style("contacthistory") != 1) error->all("Pair granular with shear history requires exactly one fix of style contacthistory");
        fix_history = static_cast<FixContactHistory*>(modify->find_fix_style("contacthistory",0));

        delete [] fixarg[3];
        delete [] fixarg;
    }
  }

  // register per-particle properties for energy tracking
  if(energytrack_enable)
  {
      if(comm->nprocs > 1) error->all("check communication settings");
      char **fixarg = new char*[9];
      if (fppaCPEn == NULL) {
        fixarg[0]=(char *) "CPEn";
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "property/atom";
        fixarg[3]=(char *) "CPEn";
        fixarg[4]=(char *) "scalar";
        fixarg[5]=(char *) "yes";
        fixarg[6]=(char *) "no";
        fixarg[7]=(char *) "no";
        fixarg[8]=(char *) "0.0";
        modify->add_fix(9,fixarg);
        fppaCPEn=static_cast<FixPropertyAtom*>(modify->find_fix_property("CPEn","property/atom","scalar",0,0));
      }
      if (fppaCDEn == NULL) {
        fixarg[0]=(char *) "CDEn";
        fixarg[1]=(char *) "all";
        fixarg[2]=(char *) "property/atom";
        fixarg[3]=(char *) "CDEn";
        fixarg[4]=(char *) "scalar";
        fixarg[5]=(char *) "yes";
        fixarg[6]=(char *) "no";
        fixarg[7]=(char *) "no";
        fixarg[8]=(char *) "0.0";
        modify->add_fix(9,fixarg);
        fppaCDEn=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEn","property/atom","scalar",0,0));
      }
      if (fppaCPEt == NULL) {
          fixarg[0]=(char *) "CPEt";
          fixarg[1]=(char *) "all";
          fixarg[2]=(char *) "property/atom";
          fixarg[3]=(char *) "CPEt";
          fixarg[4]=(char *) "scalar";
          fixarg[5]=(char *) "yes";
          fixarg[6]=(char *) "no";
          fixarg[7]=(char *) "no";
          fixarg[8]=(char *) "0.0";
          modify->add_fix(9,fixarg);
          fppaCPEt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CPEt","property/atom","scalar",0,0));
       }
       if (fppaCDEVt == NULL) {
         fixarg[0]=(char *) "CDEVt";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CDEVt";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCDEVt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEVt","property/atom","scalar",0,0));
       }
       if (fppaCDEFt == NULL) {
         fixarg[0]=(char *) "CDEFt";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CDEFt";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCDEFt=static_cast<FixPropertyAtom*>(modify->find_fix_property("CDEFt","property/atom","scalar",0,0));
       }
       if (fppaCTFW == NULL) {
         fixarg[0]=(char *) "CTFW";
         fixarg[1]=(char *) "all";
         fixarg[2]=(char *) "property/atom";
         fixarg[3]=(char *) "CTFW";
         fixarg[4]=(char *) "scalar";
         fixarg[5]=(char *) "yes";
         fixarg[6]=(char *) "no";
         fixarg[7]=(char *) "no";
         fixarg[8]=(char *) "0.0";
         modify->add_fix(9,fixarg);
         fppaCTFW=static_cast<FixPropertyAtom*>(modify->find_fix_property("CTFW","property/atom","scalar",0,0));
      }
      if (fppaDEH == NULL) {
       fixarg[0]=(char *) "DEH";
       fixarg[1]=(char *) "all";
       fixarg[2]=(char *) "property/atom";
       fixarg[3]=(char *) "DEH";
       fixarg[4]=(char *) "scalar";
       fixarg[5]=(char *) "yes";
       fixarg[6]=(char *) "no";
       fixarg[7]=(char *) "no";
       fixarg[8]=(char *) "0.0";
       modify->add_fix(9,fixarg);
       fppaDEH=static_cast<FixPropertyAtom*>(modify->find_fix_property("DEH","property/atom","scalar",0,0));
     }
    delete []fixarg;

    if(force->newton_pair == 1) error->all("Have to implement rev comm of energy terms");
  }

  // need a half neigh list and optionally a granular history neigh list

  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->half = 0;
  neighbor->requests[irequest]->gran = 1;
  if (history) {
    irequest = neighbor->request(this);
    neighbor->requests[irequest]->id = 1; 
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->granhistory = 1;
    neighbor->requests[irequest]->dnum = dnum;
  }

  // check for freeze Fix and set freeze_group_bit

  for (i = 0; i < modify->nfix; i++)
    if (strcmp(modify->fix[i]->style,"freeze") == 0) break;
  if (i < modify->nfix) freeze_group_bit = modify->fix[i]->groupbit;
  else freeze_group_bit = 0;

  // set maxrad_dynamic and maxrad_frozen for each type
  for (i = 1; i <= atom->ntypes; i++)
  onerad_dynamic[i] = onerad_frozen[i] = 0.0;

  // include future Fix pour particles as dynamic

  for (i = 0; i < modify->nfix; i++){
    for(int j=1;j<=atom->ntypes;j++)
    {
        int pour_type = 0;
        double pour_maxrad = 0.0;
        pour_type = j;
        pour_maxrad = modify->fix[i]->max_rad(pour_type);
        onerad_dynamic[pour_type] = MAX(onerad_dynamic[pour_type],pour_maxrad);
    }
  }

  //further dynamic and frozen

  double *radius = atom->radius;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (i = 0; i < nlocal; i++)
    if (mask[i] & freeze_group_bit)
      onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
    else
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,MPI_DOUBLE,MPI_MAX,world);

  init_substyle();
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void PairGran::register_compute_pair_local(ComputePairGranLocal *ptr,int &dnum_compute)
{
   if(cpl != NULL) error->all("Pair gran allows only one compute of type pair/local");
   cpl = ptr;
   dnum_compute = dnum; //history values
}

void PairGran::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   if(cpl != ptr) error->all("Illegal situation in PairGran::unregister_compute_pair_local");
   cpl = NULL;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   optional granular history list
------------------------------------------------------------------------- */

void PairGran::init_list(int id, NeighList *ptr)
{
  if (id == 0) list = ptr;
  else if (id == 1) listgranhistory = ptr;
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairGran::init_one(int i, int j)
{
  if (!allocated) allocate();

  // cutoff = sum of max I,J radii for
  // dynamic/dynamic & dynamic/frozen interactions, but not frozen/frozen

  double cutoff = maxrad_dynamic[i]+maxrad_dynamic[j];
  cutoff = MAX(cutoff,maxrad_frozen[i]+maxrad_dynamic[j]);
  cutoff = MAX(cutoff,maxrad_dynamic[i]+maxrad_frozen[j]);

  return cutoff;
}

/* ----------------------------------------------------------------------
   neighbor callback to inform pair style of neighbor list to use
   optional granular history list
------------------------------------------------------------------------- */

void PairGran::compute(int eflag, int vflag)
{
    compute(eflag,vflag,1);
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairGran::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++)
      fwrite(&setflag[i][j],sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairGran::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) fread(&setflag[i][j],sizeof(int),1,fp);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
    }
}

void PairGran::reset_dt()
{
  dt = update->dt;
}
