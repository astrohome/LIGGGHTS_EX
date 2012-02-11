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
#include "multisphere_input.h"

using namespace LAMMPS_NS;

#define MAXLINE 2048
#define DELTA 4

InputMultisphere::InputMultisphere(LAMMPS *lmp, int argc, char **argv) : Input(lmp, argc, argv)
{}

InputMultisphere::~InputMultisphere()
{}

/* ----------------------------------------------------------------------
   process clump file
------------------------------------------------------------------------- */

int InputMultisphere::clmpfile(double **xclmp,double *rclmp,int nclmps)
{
  int n;
  int iClmp = 0;

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
    if (n == 0)
      break;

    MPI_Bcast(line,n,MPI_CHAR,0,world);

    // if n = MAXLINE, line is too long
    if (n == MAXLINE) {
      char str[MAXLINE+32];
      sprintf(str,"Input line too long: %s",line);
      error->all(str);
    }

    //parse one line from the clump file
    parse_nonlammps();

    //skip empty lines
    if(narg == 0)
    {
        if (me == 0) fprintf(screen,"Note: Skipping empty line or comment line in clump file\n");
        continue;
    }

    if(iClmp >= nclmps)
        error->all("Number of clumps in file larger than number specified");

    if(narg < 4)
        error->all("Not enough arguments in one line of clump file, need to specify [xcoo ycoo zcoo radius] in each line");

    rclmp[iClmp] = atof(arg[3]);

    for(int j = 0; j < 3; j++)
       xclmp[iClmp][j] = atof(arg[j]);

    iClmp++;
  }

  return iClmp;
}

/* ----------------------------------------------------------------------
   process all input from file
------------------------------------------------------------------------- */

void InputMultisphere::clmpfile(const char *filename, double **xclmp,double *rclmp,int nclmps)
{
  if (me == 0)
  {
    nonlammps_file = fopen(filename,"r");
    if (nonlammps_file == NULL)
    {
      char str[128];
      sprintf(str,"Cannot open clump file %s",filename);
      error->one(str);
    }
  }
  else nonlammps_file = NULL;

  if(clmpfile(xclmp,rclmp,nclmps) != nclmps)
    error->all("Number of clumps in file does not match number of clumps that were specified");

  if(nonlammps_file) fclose(nonlammps_file);

}
