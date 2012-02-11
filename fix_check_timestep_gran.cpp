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

#include "string.h"
#include "stdlib.h"
#include "atom.h"
#include "update.h"
#include "math.h"
#include "error.h"
#include "fix_check_timestep_gran.h"
#include "mech_param_gran.h"
#include "fix_property_global.h"
#include "pair_gran_hooke_history.h"
#include "force.h"
#include "comm.h"
#include "modify.h"
#include "fix_wall_gran.h"
#include "fix_mesh_gran.h"
#include "stl_tri.h"
#include "neighbor.h"
#include "mympi.h"

using namespace LAMMPS_NS;

#define BIG 1000000.

/* ---------------------------------------------------------------------- */

FixCheckTimestepGran::FixCheckTimestepGran(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal fix check/timestep/gran command, not enough arguments");

  nevery = atoi(arg[3]);
  fraction_rayleigh_lim=atof(arg[4]);
  fraction_hertz_lim=atof(arg[5]);

  int iarg = 6;

  warnflag = true;
  if(iarg<narg){
      if (narg < 8) error->all("Illegal fix check/timestep/gran command, not enough arguments");
      if(strcmp(arg[iarg++],"warn")!=0) error->all("Illegal fix check/timestep/gran command, use keyword 'warn'");
      if(strcmp(arg[iarg++],"no")==0) warnflag=false;
  }

  vector_flag = 1;
  size_vector = 3;
  global_freq = nevery;
  extvector = 1;

  fraction_rayleigh = fraction_hertz = fraction_skin = 0.;
}

/* ---------------------------------------------------------------------- */

int FixCheckTimestepGran::setmask()
{
  int mask = 0;
  mask |= END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::init()
{
  //some error checks
  if(!atom->radius_flag || !atom->density_flag) error->all("Fix check/timestep/gran can only be used together with a granular atom style");

  pg = (PairGranHookeHistory*)force->pair_match("gran/hooke/history",1);
  if(!pg) pg = (PairGranHookeHistory*)force->pair_match("gran/hertz/history",1);
  if(!pg) pg = (PairGranHookeHistory*)force->pair_match("gran/hooke",1);

  if (!pg) error->all("Fix check/timestep/gran can only be used together with: gran/hooke/history, gran/hooke, gran/hertz/history");

  mpg = pg->mpg;
  int max_type = mpg->max_type();

  fwg = NULL;
  for(int ifix=0;ifix<modify->nfix;ifix++)
      if(strncmp(modify->fix[ifix]->style,"wall/gran",8)==0)
          if(static_cast<FixWallGran*>(modify->fix[ifix])->meshwall)
              fwg = static_cast<FixWallGran*>(modify->fix[ifix]);

  Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0));
  nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0));

  if(!Y || !nu) error->all("Fix check/timestep/gran only works with a pair style that defines youngsModulus and poissonsRatio");
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::end_of_step()
{
    calc_rayleigh_hertz_estims();

    double skin = neighbor->skin;
    double dt = update->dt;

    fraction_rayleigh = dt/rayleigh_time;
    fraction_hertz = dt/hertz_time;
    fraction_skin = (vmax * dt) / neighbor->skin;

    if(warnflag&&comm->me==0)
    {
        if(fraction_skin > 0.1)
        {
            if(screen)  fprintf(screen ,"WARNING: time step too large or skin too small - particles may travel a relative distance of %f per time-step, but 0.1 * skin is %f\n",vmax*dt,0.1*skin);
            if(logfile) fprintf(logfile,"WARNING: time step too large or skin too small - particles may travel a relative distance of %f per time-step, but 0.1 * skin is %f\n",vmax*dt,0.1*skin);
        }

        if(vmax * dt > r_min)
        {
            if(screen)  fprintf(screen  ,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
            if(logfile)  fprintf(logfile,"WARNING: time step way too large - particles move further than the minimum radius in one step\n");
        }

        if(fraction_rayleigh > fraction_rayleigh_lim)
        {
            if(screen) fprintf(screen,  "WARNING: time-step is %f %% of rayleigh time\n",fraction_rayleigh*100.);
            if(logfile) fprintf(logfile,"WARNING: time-step is %f %% of rayleigh time\n",fraction_rayleigh*100.);
        }
        if(fraction_hertz > fraction_hertz_lim)
        {
            if(screen) fprintf(screen,  "WARNING: time-step is %f %% of hertz time\n",fraction_hertz*100.);
            if(logfile) fprintf(logfile,"WARNING: time-step is  %f %% of hertz time\n",fraction_hertz*100.);
        }
    }
}

/* ---------------------------------------------------------------------- */

void FixCheckTimestepGran::calc_rayleigh_hertz_estims()
{
  double **v = atom->v;
  double *density = atom->density;
  double *r = atom->radius;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  int max_type = mpg->max_type();
  int min_type = 1;

  //check rayleigh time and vmax of particles
  rayleigh_time = BIG;
  r_min = BIG;
  vmax = 0;

  double vmag;
  double rayleigh_time_i;

  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit)
    {
        double shear_mod = Y->values[type[i]-1]/(2.*(nu->values[type[i]-1]+1.));
        rayleigh_time_i = M_PI*r[i]*sqrt(density[i]/shear_mod)/(0.1631*nu->values[type[i]-1]+0.8766);
        if(rayleigh_time_i < rayleigh_time) rayleigh_time = rayleigh_time_i;

        vmag = sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]+v[i][2]*v[i][2]);
        if(vmag > vmax) vmax=vmag;

        if(r[i] < r_min) r_min = r[i];
    }
  }

  MyMPI::My_MPI_Min_Scalar(r_min,world);
  MyMPI::My_MPI_Max_Scalar(vmax,world);
  MyMPI::My_MPI_Min_Scalar(rayleigh_time,world);

  //get vmax of geometry
  FixMeshGran* mesh;
  double vmax_mesh_all,vmax_mesh=0.;
  double const *vnode;
  if(fwg)
  {
      for(int imesh=0;imesh<fwg->nFixMeshGran;imesh++)
      {
          mesh=fwg->FixMeshGranList[imesh];
          if(mesh->STLdata->mesh_moving())
          {
              for(int itri=0;itri<mesh->nTri;itri++)
                  for(int inode=0;inode<3;inode++)
                  {
                      vnode=mesh->STLdata->vnode()[itri][inode];
                      vmag=sqrt(vnode[0]*vnode[0]+vnode[1]*vnode[1]+vnode[2]*vnode[2]);
                      if(vmag>vmax_mesh) vmax_mesh=vmag;
                  }
          }
      }
  }
  MPI_Allreduce(&vmax_mesh,&vmax_mesh_all,1,MPI_DOUBLE,MPI_MAX,world);
  vmax_mesh = vmax_mesh_all;

  //decide vmax - either particle-particle or particle-mesh contact
  vmax = fmax(2.*vmax,vmax+vmax_mesh);

  //check estimation for hertz time
  //this is not exact...
  //loop over all material comibinations
  //  loop all particles
  //     test collision of particle with itself
  double hertz_time_min_all,hertz_time_min = 1000000.;
  double hertz_time_i,meff,reff,Eeff;

  for (int ti= 1 ; ti < max_type+1 ; ti++)
  {
      for (int tj=  ti ; tj < max_type+1 ; tj++)
      {
          Eeff=pg->Yeff[ti][tj];

          for (int i = 0; i < nlocal; i++)
          {
            if (mask[i] & groupbit)
            {
                if(type[i]!=ti || type[i]!=tj) continue;
                meff=4.*r[i]*r[i]*r[i]*M_PI/3.*density[i];
                reff=r[i]/2.;
                hertz_time_i=2.87*pow(meff*meff/(reff*Eeff*Eeff*vmax),0.2);
                if(hertz_time_i<hertz_time_min) hertz_time_min=hertz_time_i;
            }
          }
      }
  }

  MPI_Allreduce(&hertz_time_min,&hertz_time_min_all,1,MPI_DOUBLE,MPI_MIN,world);
  hertz_time=hertz_time_min_all;
}

/* ----------------------------------------------------------------------
   return fractions of rayleigh/hertz time-step
------------------------------------------------------------------------- */

double FixCheckTimestepGran::compute_vector(int n)
{
  if(n == 0)      return fraction_rayleigh;
  else if(n == 1) return fraction_hertz;
  else if(n == 2) return fraction_skin;
  return 0.;
}
