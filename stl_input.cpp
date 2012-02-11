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

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "ctype.h"
#include "stl_input.h"
#include "style_command.h"
#include "universe.h"
#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "group.h"
#include "domain.h"
#include "output.h"
#include "thermo.h"
#include "force.h"
#include "pair.h"
#include "min.h"
#include "modify.h"
#include "compute.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "update.h"
#include "neighbor.h"
#include "special.h"
#include "variable.h"
#include "error.h"
#include "memory.h"
#include "fix_mesh_gran.h" 

using namespace LAMMPS_NS;

#define MAXLINE 2048
#define DELTA 4

InputSTL::InputSTL(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv)
{}

InputSTL::~InputSTL()
{}

/* ----------------------------------------------------------------------
   process STL file
------------------------------------------------------------------------- */

void InputSTL::stlfile(class FixMeshGran *mesh)
{
  int n;
  int *iTri = 0, *nTriMax = 0;
  int iVertex = 0;
  bool insideSolidObject = false;
  bool insideFacet = false;
  bool insideOuterLoop = false;

  iTri = &(mesh->STLdata->nTri);
  nTriMax = &(mesh->STLdata->nTriMax);

  double phix = (mesh->rot_angle[0])*M_PI/180.;
  double phiy = (mesh->rot_angle[1])*M_PI/180.;
  double phiz = (mesh->rot_angle[2])*M_PI/180.;
  double *vert_before_rot = new double[3];
  double *vert_after_rot = new double[3];

  int flag_normalize = 0;
  int flag_outside = 0;

  while (1) {
    // read one line from input script
    // if line ends in continuation char '&', concatenate next line(s)
    // n = str length of line
    if (me == 0) {
      if (fgets(line,MAXLINE,nonlammps_file) == NULL) n = 0;
      else n = strlen(line) + 1;
      while (n >= 3 && line[n-3] == '&') {
	if (fgets(&line[n-3],MAXLINE-n+3,nonlammps_file) == NULL) n = 0;
	else n = strlen(line) + 1;
      }
    }

    // bcast the line
    // if n = 0, end-of-file
    // error if label_active is set, since label wasn't encountered
    // if original input file, code is done
    // else go back to previous input file
    MPI_Bcast(&n,1,MPI_INT,0,world);
    if (n == 0) {
      break;
    }

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long
    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(str);
    }

    //parse one line from the stl file
    parse_nonlammps();

    //for (int jj=0; jj<5; jj++) fprintf(screen,"arg %d: %s  \n",jj,arg[jj]);

    //skip empty lines
    if(narg==0){
         if (me == 0) fprintf(screen,"Note: Skipping empty line in STL file\n");
      continue;
    }

    //detect begin and end of a solid object, facet and vertices
    if (strcmp(arg[0],"solid") == 0)
    {
      if (insideSolidObject) error->all("Corrupt or unknown STL file: New solid object begins without closing prior solid object.");
      insideSolidObject=true;
      if (me == 0){
         fprintf(screen,"Solid body detected in STL file\n");
       }
    }
    else if (strcmp(arg[0],"endsolid") == 0)
    {
       if (!insideSolidObject) error->all("Corrupt or unknown STL file: End of solid object found, but no begin.");
       insideSolidObject=false;
       if (me == 0) {
         fprintf(screen,"End of solid body detected in STL file.\n");
       }
    }

    //detect begin and end of a facet within a solids object
    else if (strcmp(arg[0],"facet") == 0)
    {
      if (insideFacet) error->all("Corrupt or unknown STL file: New facet begins without closing prior facet.");
      if (!insideSolidObject) error->all("Corrupt or unknown STL file: New facet begins outside solid object.");
      insideFacet=true;
      if (me == 0){
         //fprintf(screen,"  Facet detected in solid body.\n");

       }
      (*iTri)++;

      //reallocate if more than n
      if (*iTri>=*nTriMax){
        // if (me == 0)fprintf(screen,"Growing STL input arrays\n");
        mesh->STLdata->grow_arrays();
      }

      //check for keyword normal belonging to facet
      if (strcmp(arg[1],"normal")!=0) error->all("Corrupt or unknown STL file: Facet normal not defined.");

      //import facet normal, rotate it also

      for (int j=0;j<3;j++) vert_before_rot[j]=atof(arg[2+j]);
      vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
      vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
      vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

      //check if face normal is normalized
      double len=sqrt(vert_after_rot[0]*vert_after_rot[0]+vert_after_rot[1]*vert_after_rot[1]+vert_after_rot[2]*vert_after_rot[2]);
      if(fabs(1.-len)>0.000001)
      {
          flag_normalize=1;
          vert_after_rot[0]/=len;
          vert_after_rot[1]/=len;
          vert_after_rot[2]/=len;
      }
      for (int j=0;j<3;j++) mesh->STLdata->facenormal[(*iTri)-1][j]=vert_after_rot[j];
      //for (int j=0;j<3;j++) mesh->STLdata->facenormal[(*iTri)-1][j]=atof(arg[2+j]);

    }
    else if (strcmp(arg[0],"endfacet") == 0)
    {
       if (!insideFacet) error->all("Corrupt or unknown STL file: End of facet found, but no begin.");
       insideFacet=false;
       if (iVertex!=3) error->all("Corrupt or unknown STL file: Number of vertices not equal to three (no triangle).");
       if (me == 0) {
         //fprintf(screen,"  End of facet detected in in solid body.\n");
       }
    }

    //detect begin and end of an outer loop within a facet
    else if (strcmp(arg[0],"outer") == 0)
    {
      if (insideOuterLoop) error->all("Corrupt or unknown STL file: New outer loop begins without closing prior outer loop.");
      if (!insideFacet) error->all("Corrupt or unknown STL file: New outer loop begins outside facet.");
      insideOuterLoop=true;
      iVertex=0;
      if (me == 0){
         //fprintf(screen,"    Outer loop detected in facet.\n");
       }
    }
    else if (strcmp(arg[0],"endloop") == 0)
    {
       if (!insideOuterLoop) error->all("Corrupt or unknown STL file: End of outer loop found, but no begin.");
       insideOuterLoop=false;
       if (me == 0) {
         //fprintf(screen,"    End of outer loop detected in facet.\n");
       }
    }

    else if (strcmp(arg[0],"vertex") == 0)
    {
       if (!insideOuterLoop) error->all("Corrupt or unknown STL file: Vertex found outside a loop.");

       if (me == 0) {
         //fprintf(screen,"      Vertex found.\n");
       }

      double ***nd;nd=mesh->STLdata->node;

      //read the vertex, translate and scale it
      for (int j=0;j<3;j++) vert_before_rot[j]=(atof(arg[1+j])+(mesh->off_fact[j]))*(mesh->scale_fact);
      //rotate the vertex
      vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
      vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
      vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

      //store the vertex
      for (int j=0;j<3;j++) nd[(*iTri)-1][iVertex][j]=vert_after_rot[j];

      //if a vertex is outside the simulation domain, generate a warning
      if (!domain->is_in_domain(nd[(*iTri)-1][iVertex]))
      {
          flag_outside = 1;
          
      }

      iVertex++;
      if (iVertex>3) error->all("Corrupt or unknown STL file: Can not have more than 3 vertices in a facet (only triangular meshes supported).");
    }
  }

  if(flag_normalize && comm->me==0)error->warning("STL face normals were not normalized, normalizing");
  if(flag_outside && comm->me==0)error->warning("STL file may be incompatible: One or more vertices outside simulation box");

  delete []vert_before_rot;
  delete []vert_after_rot;

}

/* ----------------------------------------------------------------------
   process all input from file
------------------------------------------------------------------------- */

void InputSTL::stlfile(const char *filename, class FixMeshGran *mesh)
{

  if (me == 0) {

    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[128];
      sprintf(str,"Cannot open file %s",filename);
      error->one(str);
    }
  } else nonlammps_file = NULL;

  stlfile(mesh);

  if(nonlammps_file) fclose(nonlammps_file);

}

