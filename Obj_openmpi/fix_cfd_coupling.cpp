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
#include "modify.h"
#include "math.h"
#include "comm.h"
#include "myvector.h"
#include "fix_cfd_coupling.h"
#include "fix.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_cfd_coupling.h"
#include "cfd_regionmodel.h"
#include "style_cfd_datacoupling.h"
#include "style_cfd_regionmodel.h"

using namespace LAMMPS_NS;

#define MAXLENGTH 30

/* ---------------------------------------------------------------------- */

FixCfdCoupling::FixCfdCoupling(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  
  if(comm->me == 0 && strncmp(arg[2],"couple/cfd",10)) error->info("Some error/warning messages may appear as fix couple/cfd messages\n");

  master_flag = 1;
  for(int ii=0;ii<modify->nfix;ii++)
  {
      if(strncmp(modify->fix[ii]->style,"couple/cfd",10) == 0) master_flag = 0;
  }

  iarg = 3;

  rm = NULL;

  dc = NULL;

  if(master_flag)
  {
      
      nevery = 1;

      if (narg < 5) error->all("Illegal fix couple/cfd command");
      if(strcmp(arg[iarg],"every") && strcmp(arg[iarg],"couple_every")) error->all("Illegal fix couple/cfd command, expecting keyword 'every'");
      iarg++;

      couple_nevery = atoi(arg[iarg++]);
      if(couple_nevery < 0)error->all("Fix couple/cfd/file: every value must be >=0");

      if (0) return;
      #define CFD_DATACOUPLING_CLASS
      #define CfdDataCouplingStyle(key,Class) \
      else if (strcmp(arg[iarg],#key) == 0) dc = new Class(lmp,iarg+1,narg,arg,this);
      #include "style_cfd_datacoupling.h"
      #undef CFD_DATACOUPLING_CLASS
      else error->all("Illegal fix couple/cfd command: Unknown data coupling style - expecting 'file' or 'MPI'");

      iarg = dc->get_iarg();

      bool hasargs = true;
      while (iarg < narg && hasargs)
      {
          hasargs = false;
          if(strcmp(arg[iarg],"regionmodel") == 0)
          {
              hasargs = true;
              iarg++;
              if (0) return;
              #define CFD_REGIONMODEL_CLASS
              #define CfdRegionStyle(key,Class) \
              else if (strcmp(arg[iarg],#key) == 0) rm = new Class(lmp,iarg+1,narg,arg,this);
              #include "style_cfd_regionmodel.h"
              #undef CFD_REGIONMODEL_CLASS
              else error->all("Unknown cfd regionmodel style");
              iarg = rm->get_iarg();
          }
      }
  }

  npull = 0;
  npush = 0;
  nvalues_max = 0;

  pullnames = NULL;
  pulltypes = NULL;
  pushnames = NULL;
  pushtypes = NULL;
  pushinvoked = NULL;
  pullinvoked = NULL;
  grow_();

  ts_create = update->ntimestep;

  couple_this = 0;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::grow_()
{
      nvalues_max+=10;
      pullnames = memory->grow_2d_char_array(pullnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:valnames");
      pulltypes = memory->grow_2d_char_array(pulltypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:valtypes");
      pushinvoked = (int *) memory->srealloc(pushinvoked,MAXLENGTH*sizeof(int),"FixCfdCoupling:pushinvoked");

      pushnames = memory->grow_2d_char_array(pushnames,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushnames");
      pushtypes = memory->grow_2d_char_array(pushtypes,nvalues_max,MAXLENGTH,"FixCfdCoupling:pushtypes");
      pullinvoked = (int *) memory->srealloc(pullinvoked,MAXLENGTH*sizeof(int),"FixCfdCoupling:pullinvoked");
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

FixCfdCoupling::~FixCfdCoupling()
{
	memory->destroy_2d_char_array(pullnames);
	memory->destroy_2d_char_array(pulltypes);
	memory->destroy_2d_char_array(pushnames);
	memory->destroy_2d_char_array(pushtypes);
}

/* ---------------------------------------------------------------------- */
void FixCfdCoupling::post_create()
{
    if(dc) dc->post_create();
    special_settings();
}
/* ---------------------------------------------------------------------- */

int FixCfdCoupling::setmask()
{
  int mask = 0;
  if(master_flag && couple_nevery) mask |= END_OF_STEP;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::init()
{
  
  master = NULL;
  for(int ifix = 0; ifix < modify->nfix; ifix++)
  {
        if(strncmp("couple/cfd",modify->fix[ifix]->style,10) == 0)
        {
            if(!master)
            {
                master = static_cast<FixCfdCoupling*>(modify->fix[ifix]);
            }
        }
  }
  if(!master) master = this;

  if (strcmp(update->integrate_style,"respa") == 0)
    nlevels_respa = ((Respa *) update->integrate)->nlevels;

  if(rm) rm->init();

  init_submodel();

  if(master == this)
  {
      for(int i = 0; i < nvalues_max; i++)
      pushinvoked[i] = pullinvoked[i] = 0;
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::setup(int vflag)
{
  if (strcmp(update->integrate_style,"verlet") == 0)
    post_force(vflag);
  else {
    ((Respa *) update->integrate)->copy_flevel_f(nlevels_respa-1);
    post_force_respa(vflag,nlevels_respa-1,0);
    ((Respa *) update->integrate)->copy_f_flevel(nlevels_respa-1);
  }

  if(update->ntimestep == 0) end_of_step();

}

/* ---------------------------------------------------------------------- */
void FixCfdCoupling::push(char *name,char *type,void *&ptr)
{
    
    int found = 0;
    for(int i = 0; i < npush; i++)
    {
        if(strcmp(name,pushnames[i]) == 0 && strcmp(type,pushtypes[i]) == 0)
        {
            found = 1;
            pushinvoked[i] = 1;
        }
    }
    if(!found)
    {
        if(comm->me == 0 && screen) fprintf(screen,"LIGGGHTS could not find property %s requested by calling program. Check your model settings in LIGGGHTS.\n",name);
        lmp->error->all("This error is fatal");
    }

    return master->dc->push(name,type,ptr);
}

void FixCfdCoupling::pull(char *name,char *type,void *&ptr)
{
    
    int found = 0;
    for(int i = 0; i < npull; i++)
    {
        if(strcmp(name,pullnames[i]) == 0 && strcmp(type,pulltypes[i]) == 0)
        {
            found = 1;
            pullinvoked[i] = 1;
        }
    }
    if(!found)
    {
        if(comm->me == 0 && screen) fprintf(screen,"LIGGGHTS could not find property %s requested by calling program. Check your model settings in LIGGGHTS.\n",name);
        lmp->error->all("This error is fatal");
    }

    return master->dc->pull(name,type,ptr);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::check_datatransfer()
{
    for(int i = 0; i < npull; i++)
       if(!pullinvoked[i])
       {
            if(comm->me == 0 && screen) fprintf(screen,"Communication of property %s was not invoked by calling program, but needed by a LIGGGHTS model. Check your model settings in OF.\n",pullnames[i]);
            lmp->error->all("This error is fatal");
       }
    for(int i = 0; i < npush; i++)
       if(!pushinvoked[i])
       {
            if(comm->me == 0 && screen) fprintf(screen,"Communication of property %s was not invoked by calling program, but needed by a LIGGGHTS model. Check your model settings in OF.\n",pushnames[i]);
            lmp->error->all("This error is fatal");
       }

}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::add_pull_property(char *name,char *type)
{
   master->add_pull_prop(name,type);
}

void FixCfdCoupling::add_pull_prop(char *name,char *type)
{
    if(strlen(name) >= MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");
    if(npull >= nvalues_max) grow_();

    for(int i = 0; i < npull; i++)
    {
        if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type) == 0) return;
        if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type)) error->all("Properties added via FixCfdCoupling::add_pull_property are inconsistent");
    }

    strcpy(pullnames[npull],name);
    strcpy(pulltypes[npull],type);
    npull++;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::add_push_property(char *name,char *type)
{
   master->add_push_prop(name,type);

}

void FixCfdCoupling::add_push_prop(char *name,char *type)
{
    if(strlen(name) >= MAXLENGTH) error->all("Fix couple/cfd: Maximum string length for a variable exceeded");
    if(npush >= nvalues_max) grow_();

    for(int i = 0; i < npush; i++)
    {
        if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type) == 0) return;
        if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type)) error->all("Properties added via FixCfdCoupling::add_push_property are inconsistent");
    }

    strcpy(pushnames[npush],name);
    strcpy(pushtypes[npush],type);
    npush++;
}

/* ---------------------------------------------------------------------- */

void* FixCfdCoupling::find_pull_property(char *name,char *type,int &len1,int &len2)
{
    return find_property(0,name,type,len1,len2);
}
void* FixCfdCoupling::find_push_property(char *name,char *type,int &len1,int &len2)
{
    return find_property(1,name,type,len1,len2);
}

void* FixCfdCoupling::find_property(int push,char *name,char *type,int &len1,int &len2)
{
    
    void *ptr = NULL;
    int flag = 0;

    // check existence of property in registry
    if(push)
    {
        for(int i = 0; i < npush; i++)
            if(strcmp(pushnames[i],name) == 0 && strcmp(pushtypes[i],type) == 0) flag = 1;
    }
    else
    {
        for(int i = 0; i < npull; i++)
            if(strcmp(pullnames[i],name) == 0 && strcmp(pulltypes[i],type) == 0) flag = 1;
    }

    if(!flag) error->all("Inconsistency in FixCfdCoupling::find_property");

    // possiblility 1
    // may be atom property - look up in atom class
    if(atom->extract(name)) return atom->extract(name);

    // possiblility 2
    // may be fix property per atom - look up in modify class

    Fix *fix1 = NULL, *fix2 = NULL;

    if(strcmp(type,"scalar") == 0 || strcmp(type,"vector") == 0)
    {
       fix1 = modify->find_fix_property(name,"property/atom",type,0,0,false);
       if(fix1 && strcmp(type,"scalar") == 0) return (void*) static_cast<FixPropertyAtom*>(fix1)->vector_atom;
       if(fix1 && strcmp(type,"vector") == 0) return (void*) static_cast<FixPropertyAtom*>(fix1)->array_atom;
    }
    else if(strcmp(type,"globalscalar") == 0)
       fix2 = modify->find_fix_property(name,"property/global","scalar",0,0,false);
    else if(strcmp(type,"globalvector") == 0)
       fix2 = modify->find_fix_property(name,"property/global","vector",0,0,false);
    else if(strcmp(type,"globalmatrix") == 0)
       fix2 = modify->find_fix_property(name,"property/global","matrix",0,0,false);

    if(fix2 && strcmp(type,"globalvector") == 0)
    {
        len1 = static_cast<FixPropertyGlobal*>(fix2)->nvalues;
        return (void*) static_cast<FixPropertyGlobal*>(fix2)->values;
    }

    if(fix2 && strcmp(type,"globalarray") == 0)
    {
        len1  = static_cast<FixPropertyGlobal*>(fix2)->size_array_rows;
        len2  = static_cast<FixPropertyGlobal*>(fix2)->size_array_cols;
        return (void*) static_cast<FixPropertyGlobal*>(fix2)->array;
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::end_of_step()
{
    
    if(master != this || couple_nevery == 0) return;

    int ts = update->ntimestep;

    if((ts+1) % couple_nevery || ts_create == ts+1) couple_this = 0;
    else couple_this = 1;

    if(ts % couple_nevery || ts_create == ts) return;

    if(!dc->liggghts_is_active) return;

    if(screen && comm->me == 0) fprintf(screen,"CFD Coupling established at step %d\n",ts);

    void *dummy = NULL;

    if(ts % couple_nevery == 0)
    {
      
      if(rm) rm->rm_update();

      for(int i = 0; i < npush; i++)
      {
           
           dc->push(pushnames[i],pushtypes[i],dummy);
      }

      for(int i = 0; i < npull; i++)
      {
         
         dc->pull(pullnames[i],pulltypes[i],dummy);
      }
    }

}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixCfdCoupling::allocate_external(int    **&data, int len2,int len1,int    initvalue)
{
      master->dc->allocate_external(data,len2,len1,initvalue);
}
void FixCfdCoupling::allocate_external(double **&data, int len2,int len1,double initvalue)
{
      master->dc->allocate_external(data,len2,len1,initvalue);
}
