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
#include "group.h"
#include "comm.h"
#include "math.h"
#include "myvector.h"
#include "fix_cfd_coupling_convection.h"
#include "fix_property_atom.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::FixCfdCouplingConvection(LAMMPS *lmp, int narg, char **arg) :  FixCfdCoupling(lmp, narg, arg)
{
    fix_convectiveFlux  = fix_heatFlux = NULL;

    if(narg < iarg + 2) error->all("Fix couple/cfd/convection: Wrong number of arguments");
    if(strcmp(arg[iarg++],"T0") != 0) error->all("Fix couple/cfd/convection: Expecting keyword 'T0'");
    T0 = atof(arg[iarg++]);

    if(T0 < 0.) error->all("Fix couple/cfd/convection: T0 must be >= 0");
}

/* ---------------------------------------------------------------------- */

FixCfdCouplingConvection::~FixCfdCouplingConvection()
{

}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::pre_delete()
{
    if(fix_convectiveFlux) modify->delete_fix("convectiveHeatFlux");
}

/* ---------------------------------------------------------------------- */

int FixCfdCouplingConvection::setmask()
{
  int mask = FixCfdCoupling::setmask();
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::special_settings()
{
  
  //register convective flux
  if(!fix_convectiveFlux)
  {
        char* fixarg[11];
        fixarg[0]="convectiveHeatFlux";
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3]="convectiveHeatFlux";
        fixarg[4]="scalar"; 
        fixarg[5]="no";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fix_convectiveFlux = modify->add_fix_property_atom(9,fixarg);
  }

  //add heat transfer model if not yet active
  FixScalarTransportEquation *fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste)
  {
        char **newarg = new char*[15];
        newarg[0] = (char *) "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = (char *) "transportequation/scalar";
        newarg[3] = (char *) "equation_id";
        newarg[4] = (char *) "heattransfer";
        newarg[5] = (char *) "quantity";
        newarg[6] = (char *) "Temp";
        newarg[7] = (char *) "default_value";
        newarg[8] = new char[30];
        sprintf(newarg[8],"%f",T0);
        newarg[9] = (char *) "flux_quantity";
        newarg[10] = (char *) "heatFlux";
        newarg[11] = (char *) "source_quantity";
        newarg[12] = (char *) "heatSource";
        newarg[13] = (char *) "capacity_quantity";
        newarg[14] = (char *) "thermalCapacity";
        modify->add_fix(15,newarg);
        delete [] newarg[8];
        delete [] newarg;
  }
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::init_submodel()
{
  
  //values to send to OF
  add_push_property("Temp","scalar");

  //values to come from OF
  add_pull_property("convectiveHeatFlux","scalar");

  fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0));

  fix_convectiveFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("convectiveHeatFlux","property/atom","scalar",0,0));
}

/* ---------------------------------------------------------------------- */

void FixCfdCouplingConvection::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double *heatFlux = fix_heatFlux->vector_atom;
  double *convectiveFlux = fix_convectiveFlux->vector_atom;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      heatFlux[i] += convectiveFlux[i];
}
