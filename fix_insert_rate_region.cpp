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
#include "stdlib.h"
#include "string.h"
#include "fix_insert_rate_region.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "region.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "myvector.h"
#include "particleToInsert.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixInsertRateRegion::FixInsertRateRegion(LAMMPS *lmp, int narg, char **arg) :
  FixInsertPack(lmp, narg, arg)
{
}

/* ---------------------------------------------------------------------- */

FixInsertRateRegion::~FixInsertRateRegion()
{
}

/* ---------------------------------------------------------------------- */

void FixInsertRateRegion::calc_insertion_properties()
{
    double dt = update->dt;

    // some error checks
    if(nflowrate > 0. && massflowrate > 0.)
        error->all("Illegal fix insert/rate/region command, both 'nflowrate' and 'massflowrate' not allowed");
    if(ninsert > 0 && massinsert > 0.)
        error->all("Illegal fix insert/rate/region command, must not define both 'nparticles' and 'mass'");
    if(insert_every == 0)
        error->all("Illegal fix insert/rate/region command, 'insert_every' = once not allowed");

    // ninsert - either defined defined directly or calculated
    if(ninsert == 0)
    {
        if(massinsert > 0.) ninsert = static_cast<int>(massinsert / fix_distribution->mass_expect());
        else error->all("Illegal fix insert/pack command, must define either 'nparticles' or 'mass'");
    }

    // flow rate, ninsert_per
    if(nflowrate == 0.)
    {
        if(massflowrate == 0.) error->all("Illegal fix insert/pack command, must define either 'massrate' or 'particlerate'");
        nflowrate = massflowrate / fix_distribution->mass_expect();
    }
    else massflowrate = nflowrate * fix_distribution->mass_expect();

    ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);
}

/* ---------------------------------------------------------------------- */

int FixInsertRateRegion::calc_ninsert_this()
{
    return FixInsert::calc_ninsert_this();
}
