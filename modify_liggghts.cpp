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

#include "stdio.h"
#include "string.h"
#include "modify.h"
#include "style_compute.h"
#include "style_fix.h"
#include "atom.h"
#include "comm.h"
#include "fix.h"
#include "compute.h"
#include "group.h"
#include "update.h"
#include "domain.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   
   add a fix property
------------------------------------------------------------------------- */
FixPropertyGlobal* Modify::add_fix_property_global(int narg,char **arg)
{
    if(narg < 5) error->all("Not enough arguments to add a fix property");
    add_fix(narg,arg);
    return static_cast<FixPropertyGlobal*>(find_fix_property(arg[3],"property/global",arg[4],0,0));
}

FixPropertyAtom* Modify::add_fix_property_atom(int narg,char **arg)
{
    if(narg < 5) error->all("Not enough arguments to add a fix property");
    add_fix(narg,arg);
    return static_cast<FixPropertyAtom*>(find_fix_property(arg[3],"property/atom",arg[4],0,0));
}

/* ----------------------------------------------------------------------
   find a fix scalar transport equation
------------------------------------------------------------------------- */

FixScalarTransportEquation* Modify::find_fix_scalar_transport_equation(const char *equation_id)
{
    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,"transportequation/scalar") == 0)
      {
          FixScalarTransportEquation *fix_ste = static_cast<FixScalarTransportEquation*>(fix[ifix]);
          if(fix_ste->match_equation_id(equation_id)) return fix_ste;
      }
    return NULL;
}

/* ----------------------------------------------------------------------
   find a fix with the requested style
------------------------------------------------------------------------- */

Fix* Modify::find_fix_style_strict(const char *style, int rank)
{
    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,style) == 0)
          if(rank > 0) rank --;
          else return fix[ifix];
    return NULL;
}

Fix* Modify::find_fix_style(const char *style, int rank)
{
    int len = strlen(style);
    for(int ifix = 0; ifix < nfix; ifix++)
      if(strncmp(fix[ifix]->style,style,len) == 0)
          if(rank > 0) rank --;
          else return fix[ifix];
    return NULL;
}

/* ----------------------------------------------------------------------
   find a fix by ID, return NULL if not found
------------------------------------------------------------------------- */

Fix* Modify::find_fix_id(const char *id)
{
  int ifix;
  for (ifix = 0; ifix < nfix; ifix++)
    if (strcmp(id,fix[ifix]->id) == 0) break;
  if (ifix == nfix) return NULL;
  return fix[ifix];
}

/* ----------------------------------------------------------------------
   
   return number of fixes with requested style
------------------------------------------------------------------------- */

int Modify::n_fixes_style(const char *style)
{
    int n_fixes,len;

    n_fixes = 0;
    len = strlen(style);

    for(int ifix = 0; ifix < nfix; ifix++)
      if(strncmp(fix[ifix]->style,style,len) == 0)
          n_fixes++;

    return n_fixes;
}

int Modify::n_fixes_style_strict(const char *style)
{
    int n_fixes,len;

    n_fixes = 0;
    len = strlen(style);

    for(int ifix = 0; ifix < nfix; ifix++)
      if(strcmp(fix[ifix]->style,style) == 0)
          n_fixes++;

    return n_fixes;
}

/* ----------------------------------------------------------------------
   
   find a property registered by a fix property/global or fix property/atom
   check if it is of desired style
   return the index in the fix array
------------------------------------------------------------------------- */

Fix* Modify::find_fix_property(const char *varname,const char *style,const char *svmstyle,int len1,int len2)
{
    return find_fix_property(varname,style,svmstyle,len1,len2,true);
}

Fix* Modify::find_fix_property(const char *varname,const char *style,const char *svmstyle,int len1,int len2,bool errflag)
{
  int ifix;
  char errmsg[200];
  Fix *fix_i = NULL;

  if((strcmp(style,"property/global")) && (strcmp(style,"property/atom"))) error->all("Valid styles for find_fix_property are 'property/global' and 'property/atom'");

  if((len1 < 0) || (len2 < 0)) error->all("Lengths for find_fix_property not valid");

  for(ifix = 0; ifix < nfix; ifix++)
  {
      if(strcmp(fix[ifix]->style,"property/atom") == 0)
         fix_i = static_cast<FixPropertyAtom*>(fix[ifix])->check_fix(varname,svmstyle,len1,len2,errflag);
      else if(strcmp(fix[ifix]->style,"property/global") == 0)
         fix_i = static_cast<FixPropertyGlobal*>(fix[ifix])->check_fix(varname,svmstyle,len1,len2,errflag);

      // check_fix returns either this or NULL
      if(fix_i) return fix_i;
  }

  // no fix found
  if(errflag)
  {
      sprintf(errmsg,"Could not locate a fix/property storing value(s) for ");
      strcat(errmsg,varname);
      error->all(errmsg);
  }
  return NULL;
}

/* ----------------------------------------------------------------------
   return if fix restart in progress
------------------------------------------------------------------------- */

int Modify::fix_restart_in_progress()
{
    return  nfix_restart_global || nfix_restart_peratom;
}

