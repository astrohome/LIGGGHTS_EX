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
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "fix_insert.h"
#include "math_extra_liggghts.h"
#include "mympi.h"
#include "myvector.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001

#define LMP_DEBUGMODE_FIXINSERT false
#define LMP_DEBUG_OUT_FIXINSERT screen

/* ---------------------------------------------------------------------- */

FixInsert::FixInsert(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 7) error->all("Illegal fix insert command, not enough arguments");

  time_depend = 1; 
  restart_global = 1;

  // required args
  iarg = 3;

  if(strcmp(arg[iarg++],"seed")) error->all("Illegal fix insert command, expecting keyword 'seed'");
  seed = atoi(arg[iarg++]);
  if (seed <= 0) error->all("Illegal fix insert command, illegal seed");

  // set defaults
  init_defaults();

  xnear = NULL;

  // parse args
  
  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if(strcmp(arg[iarg],"distributiontemplate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      int ifix = modify->find_fix(arg[iarg+1]);
      if(ifix < 0 || strcmp(modify->fix[ifix]->style,"particledistribution/discrete"))
        error->all("Fix insert requires you to define a valid ID for a fix of type particledistribution/discrete");
      fix_distribution = static_cast<FixParticledistributionDiscrete*>(modify->fix[ifix]);
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"maxattempt") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      maxattempt = atoi(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"nparticles") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      ninsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"mass") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      massinsert = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"massrate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      massflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"particlerate") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      nflowrate = atof(arg[iarg+1]);
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"insert_every") == 0 || strcmp(arg[iarg],"every") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      if(strcmp(arg[iarg+1],"once") == 0) insert_every = 0;
      else insert_every = atoi(arg[iarg+1]);
      if(insert_every < 0) error->all("Illegal fix insert command, insert_every must be >= 0");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"start") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert/stream command");
      first_ins_step = atoi(arg[iarg+1]);
      if(first_ins_step < update->ntimestep + 1 && !modify->fix_restart_in_progress())
        error->all("Illegal fix insert command, 'start' step can not be before current step");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"overlapcheck") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      if(strcmp(arg[iarg+1],"yes")==0) check_ol_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) check_ol_flag = 0;
      else error->all("Illegal fix insert command");
      iarg += 2;
      hasargs = true;
    } else if (strcmp(arg[iarg],"all_in") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert command");
      if(strcmp(arg[iarg+1],"yes")==0) all_in_flag = 1;
      else if(strcmp(arg[iarg+1],"no")==0) all_in_flag = 0;
      else error->all("Illegal fix insert command");
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix insert command");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          v_insert[0] = atof(arg[iarg+2]);
          v_insert[1] = atof(arg[iarg+3]);
          v_insert[2] = atof(arg[iarg+4]);
      } else error->all("Illegal fix insert command, expecting keyword 'constant' after keyword 'vel'");
      iarg += 5;
      hasargs = true;
    }else if (strcmp(arg[iarg],"omega") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix insert command");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          omega_insert[0] = atof(arg[iarg+2]);
          omega_insert[1] = atof(arg[iarg+3]);
          omega_insert[2] = atof(arg[iarg+4]);
      } else error->all("Illegal fix insert command, expecting keyword 'constant' after keyword 'omega'");
      iarg += 5;
      hasargs = true;
    }else if (strcmp(arg[iarg],"quat") == 0) {
      if (iarg+5 > narg) error->all("Illegal fix insert command");
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          quat_insert[0] = atof(arg[iarg+2]);
          quat_insert[1] = atof(arg[iarg+3]);
          quat_insert[2] = atof(arg[iarg+4]);
          quat_insert[3] = atof(arg[iarg+5]);
      } else error->all("Illegal fix insert command, expecting keyword 'constant' after keyword 'quat'");
      iarg += 6;
      hasargs = true;
    }
    
    else if(strcmp(style,"insert") == 0) error->all("Illegal fix insert command, unknown keyword");
  }

  // default is that total # of particles to insert by this command is known
  ninsert_exists = 1;

  // default is statistically exact, non-truncated particle distributions
  truncate = 0;

  // memory not allocated initially
  ninsert_this_max = 0;

  // check for missing or contradictory settings
  sanity_check();

  //min/max type to be inserted, need that to check if material properties defined for all materials
  type_max = fix_distribution->max_type();
  type_min = fix_distribution->min_type();

  // random number generator, same for all procs
  random = new RanPark(lmp,seed);

  // allgather arrays
  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  recvcounts = new int[nprocs];
  displs = new int[nprocs];

  // set next reneighbor
  force_reneighbor = 1;
  next_reneighbor = first_ins_step;
  most_recent_ins_step = -1;

  print_stats_start_flag = 1;

  // calc max insertion radius
  int ntypes = atom->ntypes;
  maxrad = 0.;
  for(int i = 1; i <= ntypes; i++)
     maxrad = MathExtraLiggghts::max(maxrad,max_rad(i));
}

/* ---------------------------------------------------------------------- */

void FixInsert::post_create()
{
  
  // calculate ninsert, insert_every, ninsert_per
  calc_insertion_properties();

  // calc last step of insertion
  if(ninsert_exists)
  {
      if(ninsert < ninsert_per)
        final_ins_step = first_ins_step;
      else
        final_ins_step = first_ins_step + static_cast<int>(static_cast<double>(ninsert)/ninsert_per *  static_cast<double>(insert_every));
  }
  else
    final_ins_step = -1;

  // print statistics
  print_stats_start();

}

/* ---------------------------------------------------------------------- */

FixInsert::~FixInsert()
{
  delete random;
  delete [] recvcounts;
  delete [] displs;
}

/* ---------------------------------------------------------------------- */

void FixInsert::init_defaults()
{
  ninsert = ninserted = 0;
  massinsert = massinserted = 0.;
  nflowrate = massflowrate = 0.;

  insert_every = -1;
  ninsert_per = 0.;

  // 1st insertion on next timestep is default
  first_ins_step = update->ntimestep + 1;

  maxattempt = 50;

  check_ol_flag = 1;
  all_in_flag = 0;

  vectorZeroize3D(v_insert);
  vectorZeroize3D(omega_insert);

  quatUnitize4D(quat_insert);

  fix_distribution = NULL;
}

/* ---------------------------------------------------------------------- */

void FixInsert::sanity_check()
{
    if(fix_distribution == NULL) error->all("Illegal fix insert command, have to define a 'distributiontemplate'");
    if(vectorMag4DSquared(quat_insert) != 1.) error->all("Illegal fix insert command, quaternion not valid");

    if(ninsert > 0 && massinsert > 0.) error->all("Illegal fix insert command, must not define both 'nparticles' and 'mass'");
    if(nflowrate > 0. && massflowrate > 0.) error->all("Illegal fix insert command, must not define both 'particlerate' and 'massrate'");

    if(insert_every == 0 && (massflowrate > 0. || nflowrate > 0.)) error->all("Illegal fix insert/pack command, must not define 'particlerate' or 'massrate' for 'insert_every' = 0");
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_start()
{
  if (me == 0 && print_stats_start_flag && ninsert_exists) {
    if (screen)
        fprintf(screen ,"Particle insertion: %f particles every %d steps - particle rate %f  (mass rate %f)\n   %d particles (mass %f) within %d steps\n",
            ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);

    if (logfile)
        fprintf(logfile,"Particle insertion: %f particles every %d steps - particle rate %f, (mass rate %f)\n   %d particles (mass %f) within %d steps\n",
            ninsert_per,insert_every,nflowrate,massflowrate,ninsert,massinsert,final_ins_step-first_ins_step);
  }
}

/* ---------------------------------------------------------------------- */

void FixInsert::print_stats_during(int ninsert_this, double mass_inserted_this)
{
  int step = update->ntimestep;

  if (me == 0)
  {
    if (screen)
      fprintf(screen ,"Particle insertion: inserted %d particle templates (mass %f) at step %d\n - a total of %d particle templates (mass %f) inserted so far.\n",
	      ninsert_this,mass_inserted_this,step,ninserted,massinserted);

    if (logfile)
      fprintf(logfile,"Particle insertion: inserted %d particle templates (mass %f) at step %d\n - a total of %d particle templates (mass %f) inserted so far.\n",
	      ninsert_this,mass_inserted_this,step,ninserted,massinserted);
  }
}

/* ---------------------------------------------------------------------- */

int FixInsert::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsert::init()
{
    if (!atom->radius_flag || !atom->rmass_flag) error->all("Fix insert requires atom attributes radius, rmass");
    if (domain->triclinic) error->all("Cannot use fix insert with triclinic box");
    if (domain->dimension != 3) error->all("Can use fix insert for 3d simulations only");
    
    fix_rm = (FixRigidMultisphere*)(modify->find_fix_style("rigid/multisphere", 0));
}

/* ---------------------------------------------------------------------- */

int FixInsert::min_type()
{
    return type_min;
}

/* ---------------------------------------------------------------------- */

int FixInsert::max_type()
{
    return type_max;
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_rad(int type)
{
    return fix_distribution->max_rad(type);
}

/* ---------------------------------------------------------------------- */

double FixInsert::max_r_bound()
{
    return fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

double FixInsert::extend_cut_ghost()
{
    
    return 2.*fix_distribution->max_r_bound();
}

/* ---------------------------------------------------------------------- */

int FixInsert::calc_ninsert_this()
{
  if(ninsert_per == 0.) error->all("Fix insert: ninsert_per == 0.");

  // number of bodies to insert this timestep
  int ninsert_this = static_cast<int>(ninsert_per + random->uniform());
  if (ninsert_exists && ninserted + ninsert_this > ninsert) ninsert_this = ninsert - ninserted;

  return ninsert_this;
}

/* ----------------------------------------------------------------------
   perform particle insertion
------------------------------------------------------------------------- */

void FixInsert::pre_exchange()
{
  
  // just return if should not be called on this timestep
  
  if (next_reneighbor != update->ntimestep || most_recent_ins_step == update->ntimestep) return;
  most_recent_ins_step = update->ntimestep;

  // things to be done before inserting new particles
  pre_insert();

  // number of particles to insert this timestep
  int ninsert_this = calc_ninsert_this();

  // re-allocate if necessary
  if(ninsert_this > ninsert_this_max)
  {
      fix_distribution->random_init_list(ninsert_this);
      ninsert_this_max = ninsert_this;
  }

  // generate list of insertions
  // number of inserted particles can change if exact_distribution = 1
  
  ninsert_this = fix_distribution->randomize_list(ninsert_this,groupbit,truncate);
  if (ninsert_exists && ninserted + ninsert_this > ninsert) ninsert_this = ninsert - ninserted;

  if(ninsert_this == 0)
  {
      // warn if flowrate should be fulfilled
      if(nflowrate > 0. || massflowrate > 0.)
        error->warning("Particle insertion: Inserting no particle - check particle insertion settings");

      // schedule next insertion
      if (insert_every && (!ninsert_exists || ninserted < ninsert))
        next_reneighbor += insert_every;
      return;
  }
  if(ninsert_this < 0)
  {
      error->one("Particle insertion: Internal error");
  }

  // fill xnear array with particles to check overlap against
  
  // add particles in insertion volume to xnear list
  nspheres_near = 0;
  xnear = NULL;
  if(check_ol_flag)
      nspheres_near = load_xnear(ninsert_this);

  // initialize insertion counter in this step
  int ninserted_this = 0, ninserted_spheres_this = 0, ninserted_spheres_this_local = 0;
  double mass_inserted_this = 0.;

  // randomize insertion positions and set v, omega
  // also performs overlap check via xnear if requested
  // returns # bodies and # spheres that could actually be inserted
  x_v_omega(ninsert_this,ninserted_this,ninserted_spheres_this,mass_inserted_this);

  // actual particle insertion
  
  ninserted_spheres_this_local = fix_distribution->insert(ninserted_this);

  // give derived classes the chance to do some wrap-up
  finalize_insertion(ninserted_spheres_this_local);

  // give particle distributions the chance to do some wrap-up
  
  fix_distribution->finalize_insertion();

  if(ninserted_this < ninsert_this)
      error->warning("Particle insertion: Less insertions than requested");

  // set tag # of new particles beyond all previous atoms, reset global natoms
  // if global map exists, reset it now instead of waiting for comm
  // since deleting atoms messes up ghosts

  if (atom->tag_enable)
  {
    atom->tag_extend();
    atom->natoms += static_cast<double>(ninserted_spheres_this);
    if (atom->map_style)
    {
      atom->nghost = 0;
      atom->map_init();
      atom->map_set();
    }
  }

  // tally stats
  ninserted += ninserted_this;
  massinserted += mass_inserted_this;
  print_stats_during(ninserted_this,mass_inserted_this);

  // free local memory
  if(xnear) memory->destroy_2d_double_array(xnear);

  // next timestep to insert
  if (insert_every && (!ninsert_exists || ninserted < ninsert)) next_reneighbor += insert_every;
  else next_reneighbor = 0;

}

/* ----------------------------------------------------------------------
   count # of local particles that could overlap
------------------------------------------------------------------------- */

int FixInsert::count_nnear()
{
    int nlocal = atom->nlocal;
    int ncount = 0;

    for(int i = 0; i < nlocal; i++)
        ncount += overlap(i);

    return ncount;
}

/* ----------------------------------------------------------------------
   fill xnear with nearby particles
------------------------------------------------------------------------- */

int FixInsert::load_xnear(int ninsert_this)
{
  // count nearby spheres
  // setup for allgatherv
  int nspheres_near_all;
  int nspheres_near = count_nnear();

  MyMPI::My_MPI_Sum_Scalar(nspheres_near,nspheres_near_all,world);

  // data size per particle: x and radius
  int n = 4*nspheres_near;

  // preparations for allgatherv
  MPI_Allgather(&n,1,MPI_INT,recvcounts,1,MPI_INT,world);

  displs[0] = 0;
  for (int iproc = 1; iproc < nprocs; iproc++)
    displs[iproc] = displs[iproc-1] + recvcounts[iproc-1];

  // xmine is for my atoms
  // xnear is for atoms from all procs + atoms to be inserted

  double **xmine =
    memory->create_2d_double_array(nspheres_near,4,"FixInsertPack::xmine");

  xnear = memory->create_2d_double_array(nspheres_near_all + ninsert_this * fix_distribution->max_nspheres(),4,"FixInsert::xnear");

  // load up xmine array

  double **x = atom->x;
  double *radius = atom->radius;
  int nlocal = atom->nlocal;

  int ncount = 0;
  for (int i = 0; i < nlocal; i++)
    if (overlap(i)) {
      xmine[ncount][0] = x[i][0];
      xmine[ncount][1] = x[i][1];
      xmine[ncount][2] = x[i][2];
      xmine[ncount][3] = radius[i];
      ncount++;
    }

  // perform allgatherv to acquire list of nearby particles on all procs

  double *ptr = NULL;
  if (ncount) ptr = xmine[0];
  MPI_Allgatherv(ptr,4*nspheres_near,MPI_DOUBLE,xnear[0],recvcounts,displs,MPI_DOUBLE,world);

  memory->destroy_2d_double_array(xmine);

  return nspheres_near_all;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixInsert::write_restart(FILE *fp)
{
  int n = 0;
  double list[5];
  list[n++] = static_cast<double>(random->state());
  list[n++] = static_cast<double>(ninserted);
  list[n++] = static_cast<double>(first_ins_step);
  list[n++] = static_cast<double>(next_reneighbor);
  list[n++] = massinserted;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixInsert::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;
  double next_reneighbor_re;

  seed = static_cast<int> (list[n++]);
  ninserted = static_cast<int> (list[n++]);
  first_ins_step = static_cast<int> (list[n++]);
  next_reneighbor_re = static_cast<int> (list[n++]);
  massinserted = list[n++];

  random->reset(seed);

  // in order to be able to continue pouring with increased number of particles
  // if insert was already finished in run to be restarted
  if(next_reneighbor_re != 0 && ninserted < ninsert) next_reneighbor = next_reneighbor_re;
}
