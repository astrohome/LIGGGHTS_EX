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
#include "dump_localbox.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define BIG      1.0e30

/* ---------------------------------------------------------------------- */

DumpLocalBox::DumpLocalBox(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all("Illegal dump parallel/decomposition/VTK command");

  //INFO: CURRENTLY ONLY PROC 0 writes

  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multiproc!=0) error->warning("Your 'dump mesh/gran/VTK' command is writing one file per processor, where all the files contain the same data");

  format_default = NULL;

  //number of properties written out in one line with buff
  size_one=1;  //dont use buff

  lasttimestep=-1;

  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;

  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

DumpLocalBox::~DumpLocalBox()
{
  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
}

/* ---------------------------------------------------------------------- */

void DumpLocalBox::init()
{
  if (domain->triclinic == 1) error->all("Can not dump parallel/decomposition/VTK for triclinic box");
  if (binary) error->all("Can not dump parallel/decomposition/VTK in binary mode");

  // default format not needed

  delete [] format;
  format = new char[150];

  // setup function ptrs

  header_choice = &DumpLocalBox::header_item;
  pack_choice = &DumpLocalBox::pack_item;
  write_choice = &DumpLocalBox::write_item;

  // open single file, one time only

  if (multifile == 0) openfile();

  delete []xdata;
  delete []ydata;
  delete []zdata;
  delete []xdata_all;
  delete []ydata_all;
  delete []zdata_all;
  len[0] = comm->procgrid[0]+1;
  len[1] = comm->procgrid[1]+1;
  len[2] = comm->procgrid[2]+1;
  xdata = new double[len[0]];
  xdata_all = new double[len[0]];
  ydata = new double[len[1]];
  ydata_all = new double[len[1]];
  zdata = new double[len[2]];
  zdata_all = new double[len[2]];
}

/* ---------------------------------------------------------------------- */

int DumpLocalBox::modify_param(int narg, char **arg)
{
  error->warning("dump_modify keyword is not supported by 'dump parallel/decomposition/VTK' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpLocalBox::write_header(int ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpLocalBox::count()
{
  if (comm->me!=0) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

int DumpLocalBox::pack()
{
  return (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpLocalBox::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLocalBox::header_item(int ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/gran/VTK export\nASCII\n");
}

void DumpLocalBox::footer_item()
{
  return;

}

/* ---------------------------------------------------------------------- */

int DumpLocalBox::pack_item()
{
  
  xdata[0] = -BIG;
  if(comm->myloc[0] == 0) xdata[0] = domain->sublo[0];
  for(int i = 0; i < comm->procgrid[0]; i++)
  {
      xdata[i+1] = -BIG;
      if(comm->myloc[0] == i) xdata[i+1] = domain->subhi[0];
  }

  ydata[0] = -BIG;
  if(comm->myloc[1] == 0) ydata[0] = domain->sublo[1];
  for(int i = 0; i < comm->procgrid[1]; i++)
  {
      ydata[i+1] = -BIG;
      if(comm->myloc[1] == i) ydata[i+1] = domain->subhi[1];
  }

  zdata[0] = -BIG;
  if(comm->myloc[2] == 0) zdata[0] = domain->sublo[2];
  for(int i = 0; i < comm->procgrid[2]; i++)
  {
      zdata[i+1] = -BIG;
      if(comm->myloc[2] == i) zdata[i+1] = domain->subhi[2];
  }

  MPI_Allreduce(xdata,xdata_all,len[0],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(ydata,ydata_all,len[1],MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(zdata,zdata_all,len[2],MPI_DOUBLE,MPI_MAX,world);

  return 0; //nothing to pack
}

/* ---------------------------------------------------------------------- */

void DumpLocalBox::write_item(int n, double *mybuf)
{
  
  if (comm->me!=0) return;

  //ensure it is only written once in multi-proc (work-around)
  if(lasttimestep==update->ntimestep)return;
  lasttimestep=update->ntimestep;

  if(n!=0) error->all("fatal: n should be 0 in dump_mesh.cpp because pack() returns 0");

  //write the data
  fprintf(fp,"DATASET RECTILINEAR_GRID\nDIMENSIONS %d %d %d\n",len[0],len[1],len[2]);

  fprintf(fp,"X_COORDINATES %d float\n",len[0]);
  for (int i = 0; i < len[0]; i++)
     fprintf(fp,"%f ",xdata_all[i]);
  fprintf(fp,"\n");

  fprintf(fp,"Y_COORDINATES %d float\n",len[1]);
  for (int i = 0; i < len[1]; i++)
     fprintf(fp,"%f ",ydata_all[i]);
  fprintf(fp,"\n");

  fprintf(fp,"Z_COORDINATES %d float\n",len[2]);
  for (int i = 0; i < len[2]; i++)
     fprintf(fp,"%f ",zdata_all[i]);
  fprintf(fp,"\n");

  //footer not needed
}
