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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_sphere.h"
#include "atom.h"
#include "atom_vec.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "random_park.h"
#include "particleToInsert.h"
#include "region.h"
#include "domain.h"
#include "comm.h"
#include "myvector.h"
#include "fix_region_variable.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;

#define LMP_DEBUGMODE_SPHERE false

/* ---------------------------------------------------------------------- */

FixTemplateSphere::FixTemplateSphere(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (domain->dimension != 3) error->all("Fix particletemplate/sphere is for 3D simulations only");

  restart_global = 1;

  // random number generator, same for all procs
  if (narg < 4) error->all("Illegal fix particletemplate command, not enough arguments");
  seed = atoi(arg[3]);
  random = new RanPark(lmp,seed);

  iarg = 4;

  // set default values
  atom_type = 1;

  pdf_radius = NULL;
  pdf_density = NULL;

  pti = new ParticleToInsert(lmp);

  n_pti_max = 0;
  pti_list = NULL;

  reg = NULL;
  reg_var = NULL;

  //parse further args
  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"atom_type") == 0)
    {
      if (iarg+2 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
      atom_type=atoi(arg[iarg+1]);
      if (atom_type < 1) error->all("Illegal fix particletemplate command, invalid atom type (must be >=1)");
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region") == 0)
    {
      if (iarg+2 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
      int ireg = domain->find_region(arg[iarg+1]);
      if (ireg < 0) error->all("Illegal fix particletemplate command, illegal region");
      reg = domain->regions[ireg];
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"region_variable") == 0)
    {
      if (iarg+2 > narg) error->all("Illegal fix particletemplate command, not enough arguments");
      int ifix = modify->find_fix(arg[iarg+1]);
      if (ifix < 0) error->all("Illegal fix particletemplate command, illegal region/variable fix");
      reg_var = static_cast<FixRegionVariable*>(modify->fix[ifix]);
      hasargs = true;
      iarg += 2;
    }
    else if (strcmp(arg[iarg],"radius") == 0) {
      hasargs = true;
      if(strcmp(this->style,"particletemplate/sphere")) error->all("Illegal fix particletemplate command, keyword radius only valid for particletemplate/sphere");
      if (iarg+3 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
      pdf_radius = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2]);
          if( value <= 0.) error->all("Illegal fix particletemplate/sphere command, radius must be >= 0");
          pdf_radius->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
          if (iarg+4 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
          double min = atof(arg[iarg+2]);
          double max = atof(arg[iarg+3]);
          if( min <= 0. || max <= 0. || min >= max) error->all("Illegal fix particletemplate/sphere command, illegal min or max value for radius");
          pdf_radius->set_params<RANDOM_UNIFORM>(min,max);
          iarg += 4;
      }
      else if (strcmp(arg[iarg+1],"lognormal") == 0)
      {
          if (iarg+4 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if( mu <= 0. ) error->all("Illegal fix particletemplate/sphere command, illegal mu value for radius");
          if( sigma <= 0. ) error->all("Illegal fix particletemplate/sphere command, illegal sigma value for radius");
          pdf_radius->set_params<RANDOM_LOGNORMAL>(mu,sigma);
          iarg += 4;
      }
      else error->all("Illegal fix particletemplate command, invalid radius random style");
      //use mass distribution instead of number distribution
      pdf_radius->activate_mass_shift();
      volume_expect = 4.*expectancy(pdf_radius)*expectancy(pdf_radius)*expectancy(pdf_radius)*M_PI/3.;
    }
    else if (strcmp(arg[iarg],"density") == 0) {
      hasargs = true;
      if (iarg+3 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
      pdf_density = new PDF(error);
      if (strcmp(arg[iarg+1],"constant") == 0)
      {
          double value = atof(arg[iarg+2]);
          if( value <= 0.) error->all("Illegal fix particletemplate/sphere command, density must be >= 0");
          pdf_density->set_params<RANDOM_CONSTANT>(value);
          iarg += 3;
      }
      else if (strcmp(arg[iarg+1],"uniform") == 0)
      {
          if (iarg+4 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
          double min = atof(arg[iarg+2]);
          double max = atof(arg[iarg+3]);
          if( min <= 0. || max <= 0. || min >= max) error->all("Illegal fix particletemplate/sphere command, illegal min or max value for density");
          pdf_density->set_params<RANDOM_UNIFORM>(min,max);
          iarg += 4;
      }
      else if (strcmp(arg[iarg+1],"lognormal") == 0)
      {
          if (iarg+4 > narg) error->all("Illegal fix particletemplate/sphere command, not enough arguments");
          double mu = atof(arg[iarg+2]);
          double sigma = atof(arg[iarg+3]);
          if( mu <= 0. ) error->all("Illegal fix particletemplate/sphere command, illegal mu value for density");
          if( sigma <= 0. ) error->all("Illegal fix particletemplate/sphere command, illegal sigma value for density");
          pdf_density->set_params<RANDOM_LOGNORMAL>(mu,sigma);
          iarg += 4;
      }
      else error->all("Illegal fix particletemplate command, invalid density random style");
      
    }
    else if(strcmp(style,"particletemplate/sphere") == 0) error->all("Illegal fix particletemplate command, unrecognized keyword");
  }

  if(pdf_density == NULL) error->all("Illegal fix particletemplate command, have to define 'density'");

  // end here for derived classes
  if(strcmp(this->style,"particletemplate/sphere"))return;

  if(pdf_radius == NULL) error->all("Illegal fix particletemplate command, have to define 'radius'");

  //set mass and volume
  volume_expect = pow(2.*expectancy(pdf_radius),3.)*M_PI/6.;
  mass_expect = expectancy(pdf_density) * volume_expect;
}

/* ---------------------------------------------------------------------- */

FixTemplateSphere::~FixTemplateSphere()
{
    delete random;

    delete pdf_density;
    delete pdf_radius;

    if(strcmp(style,"particletemplate/sphere") == 0)
    {
        delete pti;
        if(pti_list) delete_ptilist();
    }
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::setmask()
{
  int mask = 0;
  return mask;
}

/* ----------------------------------------------------------------------*/

Region* FixTemplateSphere::region()
{
    if(reg_var) return reg_var->region();
    else return reg;
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::randomize_single()
{
    
    pti->atom_type = atom_type;

    // randomize radius
    double radius = rand(pdf_radius,random);
    pti->radius_ins[0] = pti->r_bound_ins = radius;

    // randomize density
    pti->density_ins = rand(pdf_density,random);

    // calculate volume and mass
    pti->volume_ins = radius * radius * radius * 4.*M_PI/3.;
    pti->mass_ins = pti->density_ins*pti->volume_ins;

    // init insertion position
    vectorZeroize3D(pti->x_ins[0]);

    pti->groupbit = groupbit;

}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::init_ptilist(int n_random_max)
{
    if(pti_list) error->all("invalid FixTemplateSphere::init_list()");
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");
    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsert(lmp);
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::delete_ptilist()
{
    if(n_pti_max == 0) return;

    for(int i = 0; i < n_pti_max; i++)
       delete pti_list[i];

    memory->sfree(pti_list);
    pti_list = NULL;
    n_pti_max = 0;
}

/* ----------------------------------------------------------------------*/

void FixTemplateSphere::randomize_ptilist(int n_random,int distribution_groupbit)
{
    for(int i = 0; i < n_random; i++)
    {
        
        pti_list[i]->atom_type = atom_type;

        // randomize radius
        double radius = rand(pdf_radius,random);

        pti_list[i]->radius_ins[0] = pti_list[i]->r_bound_ins = radius;

        // randomize density
        pti_list[i]->density_ins = rand(pdf_density,random);

        // calculate volume and mass
        pti_list[i]->volume_ins = radius * radius * radius * 4.*M_PI/3.;
        pti_list[i]->mass_ins = pti_list[i]->density_ins*pti_list[i]->volume_ins;

        // init insertion position
        vectorZeroize3D(pti_list[i]->x_ins[0]);
        vectorZeroize3D(pti_list[i]->v_ins);
        vectorZeroize3D(pti_list[i]->omega_ins);

        pti_list[i]->groupbit = groupbit | distribution_groupbit; 
    }
    
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::max_rad()
{
    return pdf_max(pdf_radius);
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::max_r_bound()
{
    return pdf_max(pdf_radius);
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::volexpect()
{
    if(volume_expect < 1e-12) error->all("Fix template/sphere: Volume expectancy too small");
    return volume_expect;
}

/* ----------------------------------------------------------------------*/

double FixTemplateSphere::massexpect()
{
    return mass_expect;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::number_spheres()
{
    return 1;
}

/* ----------------------------------------------------------------------*/

int FixTemplateSphere::type()
{
    return atom_type;
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixTemplateSphere::write_restart(FILE *fp)
{
  int n = 0;
  double list[1];
  list[n++] = static_cast<int>(random->state());

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixTemplateSphere::restart(char *buf)
{
  int n = 0;
  double *list = (double *) buf;

  seed = static_cast<int> (list[n++]);

  random->reset(seed);
}
