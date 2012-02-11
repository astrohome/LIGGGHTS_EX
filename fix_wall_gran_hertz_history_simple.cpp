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
#include "fix_wall_gran_hertz_history_simple.h"
#include "pair_gran_hertz_history_simple.h"
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

FixWallGranHertzHistorySimple::FixWallGranHertzHistorySimple(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistorySimple(lmp, narg, arg)
{}

/* ---------------------------------------------------------------------- */

#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void FixWallGranHertzHistorySimple::deriveContactModelParams(int ip, double deltan, double meff_wall, double &kn, double &kt, double &gamman, double &gammat, double &xmu,double &rmu)  
{
    double polyhertz = sqrt(reff_wall*deltan);

    kn = polyhertz*k_n[itype][atom_type_wall];
    kt = polyhertz*k_t[itype][atom_type_wall];
    gamman = polyhertz*meff_wall*gamma_n[itype][atom_type_wall];
    gammat = polyhertz*meff_wall*gamma_t[itype][atom_type_wall];

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
