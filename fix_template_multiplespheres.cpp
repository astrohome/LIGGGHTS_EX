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

/* -------------------------------------------------------------------------
Thanks to Chris Stoltz (P&G) for providing
a Fortran version of the MC integrator
------------------------------------------------------------------------- */

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "fix_template_multiplespheres.h"
#include "math_extra.h"
#include "math_extra_liggghts.h"
#include "myvector.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "modify.h"
#include "comm.h"
#include "output.h"
#include "memory.h"
#include "error.h"
#include "random_mars.h"
#include "random_park.h"
#include "fix_rigid.h"
#include "particleToInsert.h"
#include "multisphere_input.h"

using namespace LAMMPS_NS;
using namespace LMP_PROBABILITY_NS;

#define LARGE 1e8
#define EPSILON 1.0e-7
#define N_SHUFFLE_BOUND 200

/* ---------------------------------------------------------------------- */

FixTemplateMultiplespheres::FixTemplateMultiplespheres(LAMMPS *lmp, int narg, char **arg) :
  FixTemplateSphere(lmp, narg, arg)
{
  if(pdf_density->rand_style() != RANDOM_CONSTANT) error->all("Fix particletemplate/multiplespheres currently only supports constant density");
  if(pdf_radius) error->all("Fix particletemplate/multiplespheres currently does not support keyword 'radius'");
  if(domain->dimension != 3) error->all("Fix particletemplate/multiplespheres only supports 3D simulations");

  // parse number of spheres
  if (strcmp(arg[iarg++],"nspheres") != 0) error->all("Illegal fix particletemplate command, expecting argument 'nspheres'");
  nspheres = atoi(arg[iarg++]);
  if(nspheres < 2) error->all("Illegal fix particletemplate/multiplespheres command, illegal number of spheres");

  // allocate arrays
  x_sphere = memory->create_2d_double_array(nspheres,3,"FixTemplateMultiplespheres:x_sphere");
  r_sphere = new double[nspheres];

  // re-create pti with correct nspheres
  delete pti;
  pti = new ParticleToInsert(lmp,nspheres);

  for (int i = 0; i < 3; i++) {
      x_min[i] = LARGE;
      x_max[i] = -LARGE;
  }

  // parse args

  if (narg < iarg+4) error->all("Illegal fix particletemplate/multisphere command, not enough arguments");
  if (strcmp(arg[iarg++],"ntry") != 0) error->all("Fix particletemplate/multisphere needs ntry to be defined");
  ntry = static_cast<int>(atoi(arg[iarg++]));

  bool spheres_read = false;

  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
    hasargs = false;

    if (strcmp(arg[iarg],"spheres") == 0)
    {
      hasargs = true;
      spheres_read = true;
      iarg++;

      if (strcmp(arg[iarg],"file") == 0)
      {
          iarg++;
          if (narg < iarg+3) error->all("Illegal fix particletemplate/multisphere command, not enough arguments");

          char *clmp_filename = arg[iarg++];

          if (strcmp(arg[iarg++],"scale") != 0) error->all("Illegal fix particletemplate command, you have to specify a scale factor");
          scale_fact = atof(arg[iarg++]);
          if(scale_fact<=0.) error->all("Illegal fix particletemplate command, scale factor must be >0");

          // allocate input class, try to open file, read data from file
          InputMultisphere *myclmp_input = new InputMultisphere(lmp,0,NULL);
          myclmp_input->clmpfile(clmp_filename,x_sphere,r_sphere,nspheres);
          delete myclmp_input;

          for(int i = 0; i < nspheres; i++)
          {
              if(r_sphere[i] <= 0.) error->all("Illegal fix particletemplate command, radius must be > 0");
              r_sphere[i] *= scale_fact;
              vectorScalarMult3D(x_sphere[i],scale_fact);
          }

          // calculate bounding box
          for(int i = 0; i < nspheres; i++)
          {
              for(int j=0;j<3;j++)
              {
                if (x_sphere[i][j]-r_sphere[i]<x_min[j]) x_min[j] = x_sphere[i][j]-r_sphere[i];
                if (x_sphere[i][j]+r_sphere[i]>x_max[j]) x_max[j] = x_sphere[i][j]+r_sphere[i];
              }
          }
      }
      else
      {
          if (narg < iarg + 4*nspheres) error->all("Illegal fix particletemplate/multiplespheres command, not enough arguments");

          //read sphere r and coos, determine min and max
          for(int i = 0; i < nspheres; i++)
          {
              r_sphere[i] = atof(arg[iarg+3]);
              if(r_sphere[i] <= 0.) error->all("Illegal fix particletemplate/multiplespheres command, radius must be >0");
              for(int j = 0; j < 3; j++)
              {
                x_sphere[i][j] = atof(arg[iarg+j]);
                if (x_sphere[i][j]-r_sphere[i]<x_min[j]) x_min[j]=x_sphere[i][j]-r_sphere[i];
                if (x_sphere[i][j]+r_sphere[i]>x_max[j]) x_max[j]=x_sphere[i][j]+r_sphere[i];
              }
              iarg+=4;
          }
      }
    }
    else if(strcmp(style,"particletemplate/multiplespheres") == 0)
        error->all("Illegal fix particletemplate/multiplespheres command, unknown keyword");
  }

  if(!spheres_read) error->all("Fix particletemplate/multisphere: need to define spheres for the template");

  if(comm->me == 0 && screen) fprintf(screen,"Calculating the properties of the given template.\n   Depending on ntry, this may take a while...\n");

  if(ntry < 1e3) error->all("fix particletemplate/multisphere: ntry is too low");
  if(comm->me == 0 && ntry < 1e5) error->warning("fix particletemplate/multisphere: ntry is very low");

}

/* ---------------------------------------------------------------------- */

FixTemplateMultiplespheres::~FixTemplateMultiplespheres()
{
    memory->destroy_2d_double_array(x_sphere);
    delete []r_sphere;
}

/* ---------------------------------------------------------------------- */

void FixTemplateMultiplespheres::post_create()
{
    // calculate bounding sphere and center of mass
    // also transforms sphere coordinates so that com = 0/0/0

    calc_bounding_sphere();
    calc_center_of_mass();
}

/* ----------------------------------------------------------------------
   calc bounding sphere with iterative procedure

   do this multiple times, randomizing point order every time, see
   http://hacksoflife.blogspot.com/2009/01/randomized-bounding-spheres.html
   choose optimal result at the end - this gives linear run-time
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::calc_bounding_sphere()
{
  r_bound = LARGE;
  int *visited = new int[nspheres];
  double d[3],delta,dist;

  for(int shuffle = 0; shuffle < N_SHUFFLE_BOUND; shuffle ++)
  {
      for(int i = 0; i < nspheres; i++) visited[i] = 0;

      int isphere = -1;
      int nvisited = 0;
      double x_bound_temp[3],rbound_temp;

      while(isphere < 0 || visited[isphere] || isphere >= nspheres )
          isphere = static_cast<int>(random->uniform()*nspheres);

      nvisited++;
      visited[isphere] = 1;

      vectorCopy3D(x_sphere[isphere],x_bound_temp);
      rbound_temp = r_sphere[isphere];

      while(nvisited < nspheres)
      {
          while(isphere < 0 || visited[isphere] || isphere >= nspheres )
               isphere = static_cast<int>(random->uniform()*nspheres);

          nvisited++;
          visited[isphere] = 1;

          vectorSubtract3D(x_sphere[isphere],x_bound_temp,d);
          dist = vectorMag3D(d);

          // do nothing if sphere is completely contained in bounding sphere
          // if not contained, shift and extend bounding sphere
          if(dist + r_sphere[isphere] > rbound_temp)
          {
              double fact = (dist + r_sphere[isphere] - rbound_temp) / (2. * dist);
              vectorScalarMult3D(d,fact);
              vectorAdd3D(x_bound_temp,d,x_bound_temp);
              rbound_temp += vectorMag3D(d);
          }
          
      }
      if(rbound_temp < r_bound)
      {
          r_bound = rbound_temp;
          vectorCopy3D(x_bound_temp,x_bound);
      }
      
  }
  delete []visited;

  // do a coarse check on the validity of the bounding sphere calc
  for(int i = 0; i < nspheres; i++)
  {
      double temp[3];
      vectorSubtract3D(x_bound,x_sphere[i],temp);
      if(vectorMag3D(temp) > r_bound) error->all("Bounding sphere calculation for template failed");
  }

}

/* ----------------------------------------------------------------------
   sqr distance from x_sphere[j] to xtest
------------------------------------------------------------------------- */

double FixTemplateMultiplespheres::dist_sqr(int j, double *xtest)
{
    double dSqr = 0.;
    dSqr += (xtest[0]-x_sphere[j][0])*(xtest[0]-x_sphere[j][0]);
    dSqr += (xtest[1]-x_sphere[j][1])*(xtest[1]-x_sphere[j][1]);
    dSqr += (xtest[2]-x_sphere[j][2])*(xtest[2]-x_sphere[j][2]);
    return dSqr;
}

/* ----------------------------------------------------------------------
   generate random point in bbox
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::generate_xtry(double *x_try)
{
    x_try[0] = x_min[0]+(x_max[0]-x_min[0])*random->uniform();
    x_try[1] = x_min[1]+(x_max[1]-x_min[1])*random->uniform();
    x_try[2] = x_min[2]+(x_max[2]-x_min[2])*random->uniform();
}

/* ----------------------------------------------------------------------
   calc center of mass
------------------------------------------------------------------------- */

void FixTemplateMultiplespheres::calc_center_of_mass()
{
  // mc integration, calc volume and com, mass weight
  int nsuccess = 0;

  double x_try[3],xcm[3],dist_j_sqr;
  int n_found = 0;

  vectorZeroize3D(xcm);

  bool alreadyChecked = false;
  for(int i = 0; i < ntry; i++)
  {
      generate_xtry(x_try);

      alreadyChecked = false;
      for(int j = 0; j < nspheres; j++)
      {
          dist_j_sqr = dist_sqr(j,x_try);

          // only count once if contained in multiple spheres
          if (alreadyChecked) break;
          if(dist_j_sqr < r_sphere[j]*r_sphere[j])
          {
              xcm[0] = (xcm[0]*static_cast<double>(nsuccess)+x_try[0])/static_cast<double>(nsuccess+1);
              xcm[1] = (xcm[1]*static_cast<double>(nsuccess)+x_try[1])/static_cast<double>(nsuccess+1);
              xcm[2] = (xcm[2]*static_cast<double>(nsuccess)+x_try[2])/static_cast<double>(nsuccess+1);
              nsuccess++;
              alreadyChecked = true;
          }
      }
  }

  // expectancy values
  volume_expect = static_cast<double>(nsuccess)/static_cast<double>(ntry)*(x_max[0]-x_min[0])*(x_max[1]-x_min[1])*(x_max[2]-x_min[2]);
  mass_expect = volume_expect*expectancy(pdf_density);
  r_equiv = pow(6.*mass_expect/(8.*expectancy(pdf_density)*M_PI),1./3.);

  // transform into a system with center of mass=0/0/0

  for(int i = 0; i < nspheres; i++)
    vectorSubtract3D(x_sphere[i],xcm,x_sphere[i]);

  vectorSubtract3D(x_min,xcm,x_min);
  vectorSubtract3D(x_max,xcm,x_max);
  vectorSubtract3D(x_bound,xcm,x_bound);

}

/* ----------------------------------------------------------------------*/

double FixTemplateMultiplespheres::max_r_bound()
{
    return r_bound;
}

/* ----------------------------------------------------------------------*/

double FixTemplateMultiplespheres::max_rad()
{
    double rmax = 0.;

    for(int j = 0; j < nspheres; j++)
      if(rmax < r_sphere[j]) rmax = r_sphere[j];

    return rmax;
}

/* ----------------------------------------------------------------------*/

int FixTemplateMultiplespheres::number_spheres()
{
    return nspheres;
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::randomize_single()
{
  
  pti->nspheres = nspheres;
  pti->density_ins = expectancy(pdf_density);
  pti->volume_ins = volume_expect;
  pti->mass_ins = mass_expect;
  pti->r_bound_ins = r_bound;
  pti->atom_type = atom_type;

  for(int j = 0; j < nspheres; j++)
  {
      pti->radius_ins[j] = r_sphere[j];
      vectorCopy3D(x_sphere[j],pti->x_ins[j]);
  }

  vectorZeroize3D(pti->v_ins);
  vectorZeroize3D(pti->omega_ins);

  pti->groupbit = groupbit;
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::init_ptilist(int n_random_max)
{
    n_pti_max = n_random_max;
    pti_list = (ParticleToInsert**) memory->smalloc(n_pti_max*sizeof(ParticleToInsert*),"pti_list");

    for(int i = 0; i < n_pti_max; i++)
       pti_list[i] = new ParticleToInsert(lmp,nspheres);
}

/* ----------------------------------------------------------------------*/

void FixTemplateMultiplespheres::randomize_ptilist(int n_random,int distribution_groupbit)
{
    
    for(int i = 0; i < n_random; i++)
    {
          ParticleToInsert *pti = pti_list[i];

          pti->density_ins = expectancy(pdf_density);
          pti->volume_ins = volume_expect;
          pti->mass_ins = mass_expect;
          pti->r_bound_ins = r_bound;
          pti->atom_type = atom_type;

          for(int j = 0; j < nspheres; j++)
          {
              pti->radius_ins[j] = r_sphere[j];
              vectorCopy3D(x_sphere[j],pti->x_ins[j]);
          }

          vectorZeroize3D(pti->v_ins);
          vectorZeroize3D(pti->omega_ins);

          pti->groupbit = groupbit | distribution_groupbit; 
    }
}
