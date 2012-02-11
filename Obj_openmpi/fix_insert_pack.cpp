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
#include "fix_insert_pack.h"
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

#define SEED_OFFSET 12

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixInsertPack::FixInsertPack(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    if (strcmp(arg[iarg],"region") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert/pack command");
      int iregion = domain->find_region(arg[iarg+1]);
      if (iregion == -1) error->all("Fix insert region ID does not exist");
      ins_region = domain->regions[iregion];
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"volumefraction") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert/pack command");
      volumefraction = atof(arg[iarg+1]);
      if(volumefraction < 0. || volumefraction > 1.) error->all("Illegal fix insert/pack command, invalid volumefraction");
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"ntry_mc") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert/pack command");
      ntry_mc = atoi(arg[iarg+1]);
      if(ntry_mc < 1000) error->all("Illegal fix insert/pack command, ntry_mc must be > 1000");
      iarg += 2;
      hasargs = true;
    } else if(strcmp(style,"insert/pack") == 0) error->all("Illegal fix insert/pack command, unknown keyword");
  }

}

/* ---------------------------------------------------------------------- */

FixInsertPack::~FixInsertPack()
{

}

/* ---------------------------------------------------------------------- */

void FixInsertPack::init_defaults()
{
      volumefraction = 0.0;
      ins_region = NULL;
      ntry_mc = 100000;
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixInsertPack::calc_insertion_properties()
{
    double dt = update->dt;

    // error check on region
    if(!ins_region) error->all("Illegal fix insert/pack command, must define an insertion region");
    ins_region->reset_random(seed + SEED_OFFSET);
    region_volume = ins_region->volume_mc(ntry_mc);

    //fprintf(screen,"FixInsertPack: Volume of insertion region: %f\n",region_volume);

    if(nflowrate > 0. || massflowrate > 0.) error->all("Illegal fix insert/pack command, specifying 'nflowrate' or 'massflowrate' currently not supported");

    // ninsert
    // if ninsert not defined directly, calculate it
    if(ninsert == 0)
    {
        if(massinsert > 0. && volumefraction > 0.) error->all("Illegal fix insert command, must not define both 'mass' and 'volumefraction'");

        if(massinsert > 0.) ninsert = static_cast<int>(massinsert / fix_distribution->mass_expect());
        else if(volumefraction > 0.)
        {
            double volume_mc = ins_region->volume_mc(ntry_mc);
            ninsert = static_cast<int>(volume_mc * volumefraction / fix_distribution->vol_expect());
        }
        else error->all("Illegal fix insert/pack command, must define either 'nparticles', 'mass', or 'volumefraction'");
    }

    if(massinsert == 0.) massinsert = ninsert * fix_distribution->mass_expect();

    // insert_every and ninsert_per
    if(insert_every < 0) error->all("Illegal fix insert/pack command, must define 'insert_every'");
    if(insert_every == 0) ninsert_per = ninsert;
    if(insert_every > 0)
    {
        // flow rate, ninsert_per
        if(nflowrate == 0.)
        {
            if(massflowrate == 0.) error->all("Illegal fix insert/pack command, must define either 'massrate' or 'particlerate'");
            nflowrate = massflowrate / fix_distribution->mass_expect();
        }
        else massflowrate = nflowrate * fix_distribution->mass_expect();

        ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertPack::init()
{
  FixInsert::init();
  
}

/* ---------------------------------------------------------------------- */

inline int FixInsertPack::overlap(int i)
{
    double pos[3], rad, cut;

    vectorCopy3D(atom->x[i],pos);
    rad = atom->radius[i];

    // choose right cutoff depending on all_in_flag

    if(all_in_flag) cut = rad + maxrad;
    else cut = rad;

    if(ins_region->match_cut(pos,cut)) return 1;
    return 0;
}

/* ----------------------------------------------------------------------
   generate random positions within insertion volume
   perform overlap check via xnear if requested
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertPack::x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this)
{
    ninserted_this = ninserted_spheres_this = 0;
    mass_inserted_this = 0.;

    double pos[3];
    ParticleToInsert *pti;

    // no overlap check
    if(!check_ol_flag)
    {
        for(int itotal = 0; itotal < ninsert_this; itotal++)
        {
            pti = fix_distribution->pti_list[ninserted_this];
            double rbound = pti->r_bound_ins;

            if(all_in_flag) ins_region->generate_random_cut_away(pos,rbound);
            else ins_region->generate_random(pos);

            // could ramdonize vel, omega, quat here

            if(pos[0] == 0. && pos[1] == 0. && pos[2] == 0.) error->all("FixInsertPack::x_v_omega() illegal position");
            ninserted_spheres_this += pti->set_x_v_omega(pos,v_insert,omega_insert,quat_insert);
            mass_inserted_this += pti->mass_ins;
            ninserted_this++;

        }
    }
    // overlap check
    // account for maxattempt
    // pti checks against xnear and adds self contributions
    else
    {
        int ntry = 0;
        int maxtry = ninsert_this * maxattempt;

        while(ntry < maxtry && ninserted_this < ninsert_this)
        {
            
            pti = fix_distribution->pti_list[ninserted_this];
            double rbound = pti->r_bound_ins;

            int nins = 0;
            while(nins == 0 && ntry < maxtry)
            {
                if(all_in_flag) ins_region->generate_random_cut_away(pos,rbound);
                else ins_region->generate_random(pos);
                ntry++;

                // could ramdonize vel, omega, quat here

                nins = pti->check_near_set_x_v_omega(pos,v_insert,omega_insert,quat_insert,xnear,nspheres_near);

            }

            if(nins > 0)
            {
                ninserted_spheres_this += nins;
                mass_inserted_this += pti->mass_ins;
                ninserted_this++;
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixInsertPack::restart(char *buf)
{
    FixInsert::restart(buf);

    ins_region->reset_random(seed + SEED_OFFSET);
}
