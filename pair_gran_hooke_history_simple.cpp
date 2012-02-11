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
#include "pair_gran_hooke_history_simple.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHookeHistorySimple::PairGranHookeHistorySimple(LAMMPS *lmp) : PairGranHookeHistory(lmp)
{
    k_n = k_t = gamma_n = gamma_t = NULL;
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistorySimple::~PairGranHookeHistorySimple()
{}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistorySimple::init_substyle()
{
  int max_type = mpg->max_type();
  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties

  k_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("kn","property/global","peratomtypepair",max_type,max_type));
  k_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("kt","property/global","peratomtypepair",max_type,max_type));

  // can be either absolute damping value or relative (which is multiplied by mass afterwards)

  // first test if user specified a valid property combination
  Fix *f1,*f2,*f3,*f4;
  f1 = modify->find_fix_property("gamman","property/global","peratomtypepair",0,0,false);
  f2 = modify->find_fix_property("gammat","property/global","peratomtypepair",0,0,false);
  f3 = modify->find_fix_property("gamman_abs","property/global","peratomtypepair",0,0,false);
  f4 = modify->find_fix_property("gammat_abs","property/global","peratomtypepair",0,0,false);

  // decide for damp_massflag, check if properties have the right length
  if(f1 && f2 && !f3 && !f4)
  {
    damp_massflag = 1;
    gamma_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gamman","property/global","peratomtypepair",max_type,max_type));
    gamma_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gammat","property/global","peratomtypepair",max_type,max_type));
  }
  else if(!f1 && !f2 && f3 && f4)
  {
    damp_massflag = 0;
    gamma_n1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gamman_abs","property/global","peratomtypepair",max_type,max_type));
    gamma_t1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("gammat_abs","property/global","peratomtypepair",max_type,max_type));
  }
  else error->all("Pair gran/hooke/history/simple: Must provide either 'gamman' and 'gammat' or 'gamman_abs' and 'gammat_abs' ");

  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type));
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type));

  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          k_n[i][j] = k_n1->compute_array(i-1,j-1);
          k_t[i][j] = k_t1->compute_array(i-1,j-1);
          gamma_n[i][j] = gamma_n1->compute_array(i-1,j-1);
          gamma_t[i][j] = gamma_t1->compute_array(i-1,j-1);

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeHistorySimple::allocate_properties(int size)
{
    memory->destroy_2d_double_array(k_n);
    memory->destroy_2d_double_array(k_t);
    memory->destroy_2d_double_array(gamma_n);
    memory->destroy_2d_double_array(gamma_t);

    memory->destroy_2d_double_array(coeffFrict);
    memory->destroy_2d_double_array(coeffRollFrict);

    memory->destroy_2d_double_array(cohEnergyDens);

    k_n = memory->create_2d_double_array(size+1,size+1,"kn");
    k_t = memory->create_2d_double_array(size+1,size+1,"kt");
    gamma_n = memory->create_2d_double_array(size+1,size+1,"gamman");
    gamma_t = memory->create_2d_double_array(size+1,size+1,"gammat");

    coeffFrict = memory->create_2d_double_array(size+1,size+1,"coeffFrict");
    coeffRollFrict = memory->create_2d_double_array(size+1,size+1,"coeffRollFrict");

    cohEnergyDens = memory->create_2d_double_array(size+1,size+1,"cohEnergyDens");
}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE
inline void PairGranHookeHistorySimple::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu)
{
    kn = k_n[itype][jtype];
    kt = k_t[itype][jtype];

    if(damp_massflag)
    {
        gamman = meff*gamma_n[itype][jtype];
        gammat = meff*gamma_t[itype][jtype];
    }
    else
    {
        gamman = gamma_n[itype][jtype];
        gammat = gamma_t[itype][jtype];
    }

    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}
#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

