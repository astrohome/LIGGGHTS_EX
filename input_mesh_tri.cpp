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
#include "input.h"
#include "modify.h"
#include "update.h"
#include "error.h"
#include "domain.h"
#include "memory.h"
#include "math.h"
#include "myvector.h"
#include "input_mesh_tri.h"
#include "fix_mesh_gran.h"

using namespace LAMMPS_NS;

#define MAXLINE 2048
#define DELTA 4

InputMeshTri::InputMeshTri(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv)
{}

InputMeshTri::~InputMeshTri()
{}

/* ----------------------------------------------------------------------
   process STL file
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile(class FixMeshGran *mesh)
{
  int n;
  int iVertex = 0;

  double phix = (mesh->rot_angle[0])*M_PI/180.;
  double phiy = (mesh->rot_angle[1])*M_PI/180.;
  double phiz = (mesh->rot_angle[2])*M_PI/180.;

  int flag_normalize = 0;
  int flag_outside = 0;

  double **points;
  int ipoint,npoints = 0;
  double vert_before_rot[3], vert_after_rot[3];

  int **cells;
  int icell,ncells = 0;

  int ntris = 0;

  int iLine = 0;

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

    //parse one line from the file
    parse_nonlammps();

    //skip empty lines
    if(narg == 0){
         if (me == 0) fprintf(screen,"Note: Skipping empty line in VTK mesh file\n");
      continue;
    }

    //increase line counter
    iLine++;

    if(iLine < 3) continue;

    if(iLine == 3)
    {
        if(strcmp(arg[0],"ASCII")) error->all("Expecting ASCII VTK mesh file, cannot continue");
        continue;
    }

    if(iLine == 4)
    {
        if(strcmp(arg[0],"DATASET UNSTRUCTURED_GRID")) error->all("Expecting ASCII VTK unstructured grid mesh file, cannot continue");
        continue;
    }

    if(iLine == 5)
    {
        if(strcmp(arg[0],"POINTS")) error->all("Expecting 'POINTS' section in ASCII VTK mesh file, cannot continue");
        npoints = atoi(arg[1]);
        points = memory->create_2d_double_array(npoints,3,"input_mesh:points");
        continue;
    }

    if(iLine <= 5+npoints)
    {
        if(narg != 3) error->all("Expecting 3 values for each point in 'POINTS' section of ASCII VTK mesh file, cannot continue");

        //read the vertex, translate and scale it
        for (int j=0;j<3;j++) vert_before_rot[j]=(atof(arg[j])+(mesh->off_fact[j]))*(mesh->scale_fact);

        //rotate the vertex
        vert_after_rot[0] = vert_before_rot[0]*cos(phiy)*cos(phiz)+vert_before_rot[1]*(cos(phiz)*sin(phix)*sin(phiy)-cos(phix)*sin(phiz))+vert_before_rot[2]*(cos(phix)*cos(phiz)*sin(phiy)+sin(phix)*sin(phiz));
        vert_after_rot[1] = vert_before_rot[0]*cos(phiy)*sin(phiz)+vert_before_rot[2]*(-cos(phiz)*sin(phix)+cos(phix)*sin(phiy)*sin(phiz))+vert_before_rot[1]*(cos(phix)*cos(phiz)+sin(phix)*sin(phiy)*sin(phiz));
        vert_after_rot[2] = vert_before_rot[2]*cos(phix)*cos(phiy)+vert_before_rot[1]*cos(phiy)*sin(phix)-vert_before_rot[0]*sin(phiy);

        if (!domain->is_in_domain(vert_after_rot))
            flag_outside = 1;

        //store the vertex
        vectorCopy3D(vert_after_rot,points[ipoint]);
        ipoint++;
        continue;
    }

    if(iLine == 6+npoints)
    {
        if(strcmp(arg[0],"CELLS")) error->all("Expecting 'CELLS' section in ASCII VTK mesh file, cannot continue");
        ncells = atoi(arg[1]);
        cells = memory->create_2d_int_array(ncells,3,"input_mesh:cells");
        continue;
    }

    //copy data of all which have 3 values - can be tri, polygon etc
    if(iLine <= 6+npoints+ncells)
    {
        if(narg == 4)
            for (int j=0;j<3;j++) cells[icell][j] = atoi(arg[1+j]);
        else
            cells[icell][0] = -1;

        icell++;
        continue;
    }

    if(iLine == 7+npoints+ncells)
    {
        if(strcmp(arg[0],"CELL_TYPES")) error->all("Expecting 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
        if(ncells != atoi(arg[1]))  error->all("Inconsistency in 'CELL_TYPES' section in ASCII VTK mesh file, cannot continue");
         icell = 0;
        continue;
    }

    //only take triangles (cell type 5 according to VTK standard) - count them
    if(iLine <= 7+npoints+2*ncells)
    {
        if(strcmp(arg[0],"5")) cells[icell][0] = -1; //remove if not a tet
        else ntris++;
        icell++;
        continue;
    }

  }

  //must throw an error here since regions may not extend outside box
  if(flag_outside)error->all("VTK mesh file is incompatible with simulation box: One or more vertices outside simulation box");

  //now that everything is parsed, write the data into the mesh
  while(mesh->STLdata->nTriMax < ntris) mesh->STLdata->grow_arrays();

  int iTri;
  for(int i = 0; i < ncells; i++)
  {
      iTri = mesh->nTri;
      if(cells[i][0] == -1) continue;
      for (int j = 0; j < 3; j++)
        vectorCopy3D(points[cells[i][j]],mesh->STLdata->node[iTri][j]);
      error->all("must set normal vector");
      mesh->nTri++;
  }

  memory->destroy_2d_double_array(points);
  memory->destroy_2d_int_array(cells);
}

/* ----------------------------------------------------------------------
   process all input from filename
------------------------------------------------------------------------- */

void InputMeshTri::meshtrifile(const char *filename, class FixMeshGran *mesh)
{

  // error if another nested file still open
  // if single open file is not stdin, close it
  // open new filename and set stl___file

  if (me == 0) {

    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL) {
      char str[128];
      sprintf(str,"Cannot open mesh file %s",filename);
      error->one(str);
    }
  } else nonlammps_file = NULL;

  meshtrifile(mesh);

  if(nonlammps_file) fclose(nonlammps_file);

}

