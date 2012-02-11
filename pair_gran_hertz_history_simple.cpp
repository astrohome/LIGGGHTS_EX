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
#include "pair_gran_hertz_history_simple.h"
#include "atom.h"
#include "force.h"
#include "memory.h"
#include "modify.h"
#include "error.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairGranHertzHistorySimple::PairGranHertzHistorySimple(LAMMPS *lmp) : PairGranHookeHistorySimple(lmp)
{}

/* ---------------------------------------------------------------------- */

PairGranHertzHistorySimple::~PairGranHertzHistorySimple()
{}

/* ----------------------------------------------------------------------
 return appropriate params
------------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE
inline void PairGranHertzHistorySimple::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu)
{
    double polyhertz = sqrt(ri*rj/(ri+rj)*deltan);
    kn = polyhertz*k_n[itype][jtype];
    kt = polyhertz*k_t[itype][jtype];
    gamman = polyhertz*meff*gamma_n[itype][jtype];
    gammat = polyhertz*meff*gamma_t[itype][jtype];

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
