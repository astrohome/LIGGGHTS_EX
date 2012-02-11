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
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran_hooke_history_simple.h"
#include "pair_gran_hooke_history_simple.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixWallGranHookeHistorySimple::FixWallGranHookeHistorySimple(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistory(lmp, narg, arg)
{
    k_n1 = k_t1 = gamma_n1 = gamma_t1 = NULL;
}

/* ---------------------------------------------------------------------- */

void FixWallGranHookeHistorySimple::init_substyle()
{
  k_n = ((PairGranHookeHistorySimple*)pairgran)->k_n;
  k_t = ((PairGranHookeHistorySimple*)pairgran)->k_t;
  gamma_n = ((PairGranHookeHistorySimple*)pairgran)->gamma_n;
  gamma_t = ((PairGranHookeHistorySimple*)pairgran)->gamma_t;
  coeffFrict = ((PairGranHookeHistorySimple*)pairgran)->coeffFrict;
  coeffRollFrict = ((PairGranHookeHistorySimple*)pairgran)->coeffRollFrict;

  damp_massflag = ((PairGranHookeHistorySimple*)pairgran)->damp_massflag;

  //need to check properties for rolling friction and cohesion energy density here
  //since these models may not be active in the pair style
  
  FixPropertyGlobal *coeffRollFrict1, *cohEnergyDens1;
  int max_type = pairgran->mpg->max_type();
  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type));
  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);
          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
      }
  }
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void FixWallGranHookeHistorySimple::deriveContactModelParams(int ip, double deltan,double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu)  
{

    kn = k_n[itype][atom_type_wall];
    kt = k_t[itype][atom_type_wall];

    if(damp_massflag)
    {
        gamman = meff_wall*gamma_n[itype][atom_type_wall];
        gammat = meff_wall*gamma_t[itype][atom_type_wall];
    }
    else
    {
        gamman = gamma_n[itype][atom_type_wall];
        gammat = gamma_t[itype][atom_type_wall];
    }

    xmu=coeffFrict[itype][atom_type_wall];
    if(rollingflag)rmu=coeffRollFrict[itype][atom_type_wall];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE
