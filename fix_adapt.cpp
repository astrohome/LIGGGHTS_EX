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
#include "fix_adapt.h"
#include "atom.h"
#include "force.h"
#include "pair.h"
#include "input.h"
#include "variable.h"
#include "error.h"
#include "fix_property_atom.h"
#include "modify.h"
#include "update.h"

using namespace LAMMPS_NS;

enum{PAIR,ATOM};
enum{DIAMETER};

/* ---------------------------------------------------------------------- */

FixAdapt::FixAdapt(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 4) error->all("Illegal fix adapt command");
  nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix adapt command");

  // count # of adaptations

  nadapt = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all("Illegal fix adapt command");
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix adapt command");
      nadapt++;
      iarg += 3;
    } else error->all("Illegal fix adapt command");
  }

  // allocate per-adapt vectors

  which = new int[nadapt];
  pair = new char*[nadapt];
  param = new char*[nadapt];
  ilo = new int[nadapt];
  ihi = new int[nadapt];
  jlo = new int[nadapt];
  jhi = new int[nadapt];
  var = new char*[nadapt];
  ivar = new int[nadapt];
  pairptr = new Pair*[nadapt];
  pairindex = new int[nadapt];
  awhich = new int[nadapt];

  // parse keywords

  nadapt = 0;

  iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"pair") == 0) {
      if (iarg+6 > narg) error->all("Illegal fix adapt command");
      which[nadapt] = PAIR;
      int n = strlen(arg[iarg+1]) + 1;
      pair[nadapt] = new char[n];
      strcpy(pair[nadapt],arg[iarg+1]);
      n = strlen(arg[iarg+2]) + 1;
      param[nadapt] = new char[n];
      strcpy(param[nadapt],arg[iarg+2]);
      force->bounds(arg[iarg+3],atom->ntypes,ilo[nadapt],ihi[nadapt]);
      force->bounds(arg[iarg+4],atom->ntypes,jlo[nadapt],jhi[nadapt]);
      n = strlen(arg[iarg+5]) + 1;
      var[nadapt] = new char[n];
      strcpy(var[nadapt],arg[iarg+5]);
      nadapt++;
      iarg += 6;
    } else if (strcmp(arg[iarg],"atom") == 0) {
      if (iarg+3 > narg) error->all("Illegal fix adapt command");
      which[nadapt] = ATOM;
      int n = strlen(arg[iarg+1]) + 1;
      param[nadapt] = new char[n];
      strcpy(param[nadapt],arg[iarg+1]);
      
      if (strcmp(param[nadapt],"diameter") == 0)	rad_mass_vary_flag = 1;
      n = strlen(arg[iarg+2]) + 1;
      var[nadapt] = new char[n];
      strcpy(var[nadapt],arg[iarg+2]);
      nadapt++;
      iarg += 3;
    } else error->all("Illegal fix adapt command");
  }

  fppat=NULL; 

  if(rad_mass_vary_flag) force_reneighbor = 1;

}

/* ---------------------------------------------------------------------- */

FixAdapt::~FixAdapt()
{
  for (int i = 0; i < nadapt; i++) {
    if (which[i] == PAIR) delete [] pair[i];
    delete [] param[i];
    delete [] var[i];
  }
  delete [] which;
  delete [] pair;
  delete [] param;
  delete [] ilo;
  delete [] ihi;
  delete [] jlo;
  delete [] jhi;
  delete [] var;
  delete [] ivar;
  delete [] pairptr;
  delete [] pairindex;
  delete [] awhich;
}

/* ---------------------------------------------------------------------- */

void FixAdapt::pre_delete()
{
    if (fppat) modify->delete_fix(fixid); 
}

/* ---------------------------------------------------------------------- */

int FixAdapt::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixAdapt::post_create()
{
  // error checks

  atomstyle_flag = 0;

  for (int m = 0; m < nadapt; m++) {
    if (which[m] == PAIR) {
      pairptr[m] = force->pair_match(pair[m],1);
      if (pairptr[m] == NULL)
	error->all("Fix adapt pair style does not exist");
      pairindex[m] =
	pairptr[m]->pre_adapt(param[m],ilo[m],ihi[m],jlo[m],jhi[m]);
      if (pairindex[m] == -1)
	error->all("Fix adapt pair parameter is not recognized");
      if (pairindex[m] == -2)
	error->all("Fix adapt pair types are not valid");

    } else if (which[m] == ATOM) {
      if (strcmp(param[m],"diameter") == 0) {
	awhich[m] = DIAMETER;
	if (!atom->radius_flag)
	  error->all("Fix adapt requires atom attribute radius");
      } else error->all("Fix adapt atom attribute is not recognized");
    }

    ivar[m] = input->variable->find(var[m]);
    if (ivar[m] < 0) error->all("Variable name for fix adapt does not exist");
    if (!input->variable->equalstyle(ivar[m])&&!input->variable->atomstyle(ivar[m]))
      error->all("Variable for fix adapt is not of style equal or atom");

    if(input->variable->atomstyle(ivar[m])) atomstyle_flag = 1;

  }

  if (fppat == NULL && atomstyle_flag)
  {
    char **fixarg;
    fixarg=new char*[9];
    for (int kk=0;kk<9;kk++) fixarg[kk]=new char[50];
    
    sprintf(fixarg[0],"adaptProp_%s",id);
    sprintf(fixid,    "adaptProp_%s",id);
    fixarg[1]="all";
    fixarg[2]="property/atom";
    sprintf(fixarg[3],"adaptProp_%s",id);
    fixarg[4]="scalar"; 
    fixarg[5]="no";    
    fixarg[6]="no";    
    fixarg[7]="no";    
    fixarg[8]="0.";     
    modify->add_fix(9,fixarg);
    delete []fixarg;

  }

  if(rad_mass_vary_flag) next_reneighbor = update->ntimestep + 1;
}

void FixAdapt::init()
{

  // set params to values for initial force calculation
  // needs to happen here in init() instead of setup()
  // because modify->setup() is called after pre-Verlet forces are computed

  pre_force(0);
}

/* ---------------------------------------------------------------------- */

void FixAdapt::pre_force(int vflag)
{
  if (update->ntimestep % nevery) return;

  modify->clearstep_compute();

  if(atomstyle_flag) fppat = static_cast<FixPropertyAtom*>(modify->find_fix_property(fixid,"property/atom","scalar",0,0));

  for (int m = 0; m < nadapt; m++) {
    double value;

    if (input->variable->equalstyle(ivar[m])) value = input->variable->compute_equal(ivar[m]);
    else if(input->variable->atomstyle(ivar[m])) input->variable->compute_atom(ivar[m], igroup, fppat->vector_atom, 1, 0);
    else error->all("Wrong variable style in fix adapt - must use a scalar property");

    bool is_atomstyle = (input->variable->atomstyle(ivar[m])) == 1;

    if (which[m] == PAIR)
      pairptr[m]->adapt(pairindex[m],ilo[m],ihi[m],jlo[m],jhi[m],value);

    else if (which[m] == ATOM) {

      // set radius from diameter
      // set rmass if both rmass and density are defined

      if (awhich[m] == DIAMETER) {
        int mflag = 0;
        if (atom->rmass_flag && atom->density_flag) mflag = 1;
        double PI = 4.0*atan(1.0);

        double *radius = atom->radius;
        double *rmass = atom->rmass;
        double *density = atom->density;
        int *mask = atom->mask;
        int nlocal = atom->nlocal;

        for (int i = 0; i < nlocal; i++)
          if (mask[i] & groupbit) {
            if (is_atomstyle) value = fppat->vector_atom[i];
            
            radius[i] = 0.5*value;
            rmass[i] = 4.0*PI/3.0 * radius[i]*radius[i]*radius[i] * density[i];
          }
        
      }
    }
  }
  modify->addstep_compute(update->ntimestep + nevery);
}
