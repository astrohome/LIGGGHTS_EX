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
#include "lammps.h"
#include "atom.h"
#include "mpi.h"
#include "math.h"
#include "modify.h"
#include "mech_param_gran.h"
#include "error.h"
#include "fix_property_global.h"
#include "memory.h"
#include "fix_insert.h"
#include "fix_pour.h"
#include "fix_pour_dev.h"
#include "fix_wall_gran_hooke_history.h"
#include "fix_mesh_gran.h"

using namespace LAMMPS_NS;

MechParamGran::MechParamGran(LAMMPS *lmp): Pointers(lmp)
{
}

MechParamGran::~MechParamGran()
{
}

int MechParamGran::max_type()
{
  //loop over all particles to check how many atom types are present
  mintype=1;
  maxtype=1;

  for (int i=0;i<atom->nlocal;i++)
  {
      if (atom->type[i]<mintype) mintype=atom->type[i];
      if (atom->type[i]>maxtype) maxtype=atom->type[i];
  }

  //check all fixes of type pour
  for(int i=0;i<modify->nfix;i++)
  {
      
      if(strncmp(modify->fix[i]->style,"insert",6)==0)
      {
          //fprintf(
          int tp_min=static_cast<FixInsert*>(lmp->modify->fix[i])->min_type();
          int tp_max=static_cast<FixInsert*>(lmp->modify->fix[i])->max_type();
          if(tp_min<mintype) mintype=tp_min;
          if(tp_max>maxtype) maxtype=tp_max;
      }
      else if(strncmp(modify->fix[i]->style,"pour/dev",7)==0||strcmp(lmp->modify->fix[i]->style,"pour/multisphere")==0)
      {
          int tp_min=static_cast<FixPourDev*>(lmp->modify->fix[i])->ntype_min;
          int tp_max=static_cast<FixPourDev*>(lmp->modify->fix[i])->ntype_max;
          if(tp_min<mintype) mintype=tp_min;
          if(tp_max>maxtype) maxtype=tp_max;
      }
      else if(strncmp(modify->fix[i]->style,"pour",4)==0)
      {
          int tp=static_cast<FixPour*>(lmp->modify->fix[i])->ntype;
          if(tp<mintype) mintype=tp;
          if(tp>maxtype) maxtype=tp;
      }
      else if(strncmp(modify->fix[i]->style,"wall/gran",8)==0)
      {
          FixWallGran* fwg=static_cast<FixWallGran*>(lmp->modify->fix[i]);
          if(fwg->meshwall)
          {
              for(int j=0;j<fwg->nFixMeshGran;j++)
              {
                int tp=fwg->FixMeshGranList[j]->atom_type_wall;
                if(tp<mintype) mintype=tp;
                if(tp>maxtype) maxtype=tp;
              }
          }
          else
          {
            int tp=fwg->atom_type_wall;
            if(tp<mintype) mintype=tp;
            if(tp>maxtype) maxtype=tp;
          }
      }
  }

  //Get min/max from other procs
  int mintype_all,maxtype_all;
  MPI_Allreduce(&mintype,&mintype_all, 1, MPI_INT, MPI_MIN, world);
  MPI_Allreduce(&maxtype,&maxtype_all, 1, MPI_INT, MPI_MAX, world);
  mintype=mintype_all;
  maxtype=maxtype_all;

  //error check
  if(mintype != 1) error->all("Atom types must start from 1 for granular simulations");
  if(maxtype > atom->ntypes) error->all("Please increase the number of atom types in the 'create_box' command to match the number of atom types you use in the simulation");

  return maxtype;
}
