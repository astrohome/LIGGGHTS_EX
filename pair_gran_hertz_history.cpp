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
#include "stdio.h"
#include "string.h"
#include "pair_gran_hertz_history.h"
#include "atom.h"
#include "force.h"
#include "fix_property_global.h"
#include "error.h"
#include "modify.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHertzHistory::PairGranHertzHistory(LAMMPS *lmp) :
  PairGranHookeHistory(lmp)
{
    charVelflag = 0;
}

/* ----------------------------------------------------------------------
   contact model parameters derived for hertz model 
------------------------------------------------------------------------- */

inline void PairGranHertzHistory::deriveContactModelParams(int &ip, int &jp,double &meff, double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu) 
{
    #define LMP_GRAN_DEFS_DEFINE
    #include "pair_gran_defs.h"
    #undef LMP_GRAN_DEFS_DEFINE

    double reff=ri*rj/(ri+rj);

    double sqrtval = sqrt(reff*deltan);

    double Sn=2.*Yeff[itype][jtype]*sqrtval;
    double St=8.*Geff[itype][jtype]*sqrtval;

    kn=4./3.*Yeff[itype][jtype]*sqrtval;
    kt=St;
    gamman=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(Sn*meff);
    gammat=-2.*sqrtFiveOverSix*betaeff[itype][jtype]*sqrt(St*meff);
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];

    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;
    /* NP
    char *testmsg=new char[200];
    sprintf(testmsg,"Yeff=%f,reff=%f,deltan=%f, kn=%f, kt=%f, gamman=%f, gammat=%f, xmu=%f\n",Yeff, reff,deltan, kn,kt,gamman,gammat,xmu);
    error->warning(testmsg);
    delete []testmsg;*/

    #define LMP_GRAN_DEFS_UNDEFINE
    #include "pair_gran_defs.h"
    #undef LMP_GRAN_DEFS_UNDEFINE

    return;
}

