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
#include "dump_mesh.h"
#include "domain.h"
#include "atom.h"
#include "update.h"
#include "group.h"
#include "error.h"
#include "fix.h"
#include "fix_mesh_gran.h"
#include "modify.h"
#include "comm.h"

using namespace LAMMPS_NS;

#define nFixMeshGran 20 //a maximum of nFixMeshGran fixes of type 'mesh/gran' is allowed

/* ---------------------------------------------------------------------- */

DumpMesh::DumpMesh(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg < 6) error->all("Illegal dump mesh/gran/VTK command");

  //info: CURRENTLY ONLY PROC 0 writes

  //multifile=1;             // 0 = one big file, 1 = one file per timestep
  //multiproc=0;             // 0 = proc 0 writes for all, 1 = one file/proc
  if (multiproc) error->warning("Your 'dump mesh/gran/VTK' command is writing one file per processor, where all the files contain the same data");
  if (!multifile) error->all("Illegal dump mesh/gran/VTK command, you have to write one file per time-step by using the '*' in the file-name");

  STLList=new STLtri*[nFixMeshGran];
  for (int i=0;i<nFixMeshGran;i++) STLList[i]=(STLtri*)NULL;

  format_default = NULL;

  //number of properties written out in one line with buff
  size_one=1;  //dont use buff

  lasttimestep=-1;

  int iarg = 5;
  dump_what = 0;

  bool hasargs = true;
  while (iarg < narg && hasargs)
  {
      hasargs = false;
      if(strcmp(arg[iarg],"stress")==0)
      {
          dump_what |= DUMP_STRESS;
          iarg++;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"stresscomponents")==0)
      {
          dump_what |= DUMP_STRESSCOMPONENTS;
          iarg++;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"id")==0)
      {
          dump_what |= DUMP_ID;
          iarg++;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"wear")==0)
      {
          dump_what |= DUMP_WEAR;
          iarg++;
          hasargs = true;
      }
      else if(strcmp(arg[iarg],"vel")==0)
      {
          dump_what |= DUMP_VEL;
          iarg++;
          hasargs = true;
      }
  }

  if(dump_what == 0) error->all("Dump mesh: No dump quantity selected");

  meshid_len=narg-iarg;
  meshid=new char*[meshid_len];
  for (int k=0;k<meshid_len;k++) meshid[k]=new char[40];

  int k=0;
  while(iarg<narg)
      strcpy(meshid[k++],arg[iarg++]);
}

/* ---------------------------------------------------------------------- */

DumpMesh::~DumpMesh()
{
    for (int k=0;k<meshid_len;k++) delete []meshid[k];
}

/* ---------------------------------------------------------------------- */

void DumpMesh::init()
{
  if (domain->triclinic == 1) error->all("Can not dump STL files for triclinic box");
  if (binary) error->all("Can not dump STL files in binary mode");

  STLList_len=0;
  for (int ifix = 0; ifix < modify->nfix; ifix++)
  {
    if (strncmp(modify->fix[ifix]->style,"mesh/gran",8) == 0)  {

        bool found = false;
        for (int k=0;k<meshid_len;k++)
            if(strcmp(modify->fix[ifix]->id,meshid[k])==0) found=true;
        if(!found && meshid_len > 0) continue;

        STLList[STLList_len]=static_cast<FixMeshGran*>(modify->fix[ifix])->STLdata;
        if(dump_what==DUMP_STRESS && !static_cast<FixMeshGran*>(modify->fix[ifix])->analyseStress) error->warning("You want to dump stress, but the stress is not calculated - use fix mesh/gran/stressanalysis");
        STLList_len++;
    }
    if (STLList_len==19) error->all("Dump mesh/gran/VTK can not process all fixes of type mesh/gran, boost nFixMeshGran in dump_mesh.cpp and recompile");
  }

  if(STLList_len<meshid_len) error->warning("Dump stl could not locate all fixes of type 'mesh/gran' that you wish to dump");
  if (STLList_len==0) error->warning("Dump stl cannot find any fix of type 'mesh/gran' to dump");

  for (int i = 0; i < STLList_len; i++){
      if(!STLList[i]->mesh_moving()) error->warning("Dump mesh/gran/VTK contains a mesh with no velocity allocated, will not dump velocity");
  }

  // default format not needed

  delete [] format;
  format = new char[150];

  // setup function ptrs

  header_choice = &DumpMesh::header_item;
  pack_choice = &DumpMesh::pack_item;
  write_choice = &DumpMesh::write_item;

  // open single file, one time only

  if (multifile == 0) openfile();
}

/* ---------------------------------------------------------------------- */

int DumpMesh::modify_param(int narg, char **arg)
{
  error->warning("dump_modify keyword is not supported by 'dump mesh/gran/VTK' and is thus ignored");
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpMesh::write_header(int ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

int DumpMesh::count()
{
  if (comm->me!=0) return 0;
  return 1;
}

/* ---------------------------------------------------------------------- */

int DumpMesh::pack()
{
  return (this->*pack_choice)();
}

/* ---------------------------------------------------------------------- */

void DumpMesh::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpMesh::header_item(int ndump)
{
  if (comm->me!=0) return;
  fprintf(fp,"# vtk DataFile Version 2.0\nLIGGGHTS mesh/gran/VTK export\nASCII\n");
}

void DumpMesh::footer_item()
{
  return;
  
}

/* ---------------------------------------------------------------------- */

int DumpMesh::pack_item()
{
  return 0; //nothing to pack
}

/* ---------------------------------------------------------------------- */

void DumpMesh::write_item(int n, double *mybuf) 
{
  if (comm->me!=0) return;

  //ensure it is only written once in multi-proc (work-around)
  if(lasttimestep==update->ntimestep)return;
  lasttimestep=update->ntimestep;

  if(n!=0) error->all("fatal: n should be 0 in dump_mesh.cpp because pack() returns 0");

  int nTriGes=0;
  for (int i = 0; i < STLList_len; i++){
      nTriGes+=STLList[i]->nTri;
  }
  //write point data
  fprintf(fp,"DATASET POLYDATA\nPOINTS %d float\n",3*nTriGes);
  double x,y,z;
  for (int i = 0; i < STLList_len; i++) {
      for (int j=0;j<STLList[i]->nTri;j++){
          for(int k=0;k<3;k++){
              x=STLList[i]->node[j][k][0];
              y=STLList[i]->node[j][k][1];
              z=STLList[i]->node[j][k][2];
              fprintf(fp,"%f %f %f\n",x,y,z);
          }
      }
  }
  //write polygon data
  fprintf(fp,"POLYGONS %d %d\n",nTriGes,4*nTriGes);
  int m=0;
  for (int i = 0; i < nTriGes; i++) {
      fprintf(fp,"%d %d %d %d\n",3,m,m+1,m+2);
      m+=3;
  }

  fprintf(fp,"CELL_DATA %d\n",nTriGes);

  if(dump_what & DUMP_STRESS)
  {
      //write pressure and shear stress
      fprintf(fp,"SCALARS pressure float 1\nLOOKUP_TABLE default\n");
      for (int i = 0; i < STLList_len; i++) {
          for (int j=0;j<STLList[i]->nTri;j++){
              fprintf(fp,"%f\n",fabs((STLList[i]->fn_fshear[j][0])/(STLList[i]->Area[j])));
              //fprintf(fp,"%f\n",(STLList[i]->Area[j]));
          }
      }
      fprintf(fp,"SCALARS shearstress float 1\nLOOKUP_TABLE default\n");
      for (int i = 0; i < STLList_len; i++) {
          for (int j=0;j<STLList[i]->nTri;j++){
              fprintf(fp,"%f\n",fabs((STLList[i]->fn_fshear[j][1])/(STLList[i]->Area[j])));
              //fprintf(fp,"%f\n",(STLList[i]->Area[j]));
          }
      }
  }
  if(dump_what & DUMP_ID)
  {
      //write id
      fprintf(fp,"SCALARS meshid int 1\nLOOKUP_TABLE default\n");
      for (int i = 0; i < STLList_len; i++) {
          for (int j=0;j<STLList[i]->nTri;j++){
              fprintf(fp,"%d\n",j);
              //fprintf(fp,"%f\n",(STLList[i]->Area[j]));
          }
      }
  }
  if(dump_what & DUMP_WEAR)
  {
      //write wear data
      fprintf(fp,"SCALARS wear float 1\nLOOKUP_TABLE default\n");
      for (int i = 0; i < STLList_len; i++) {
          for (int j=0;j<STLList[i]->nTri;j++){
              fprintf(fp,"%f\n",STLList[i]->wear[j]);
              //fprintf(fp,"%f\n",(STLList[i]->Area[j]));
          }
      }
  }

  if(dump_what & DUMP_VEL)
  {
      //write vel data
      fprintf(fp,"VECTORS v float 3\nLOOKUP_TABLE default\n");
      for (int i = 0; i < STLList_len; i++) {
          int vel_flag = STLList[i]->mesh_moving();
          for (int j=0;j<STLList[i]->nTri;j++){
              if(vel_flag)
              {
                  double vx =  (STLList[i]->v_node[j][0][0] + STLList[i]->v_node[j][1][0] + STLList[i]->v_node[j][2][0]) / 3.;
                  double vy =  (STLList[i]->v_node[j][0][1] + STLList[i]->v_node[j][1][1] + STLList[i]->v_node[j][2][1]) / 3.;
                  double vz =  (STLList[i]->v_node[j][0][2] + STLList[i]->v_node[j][1][2] + STLList[i]->v_node[j][2][2]) / 3.;
                  fprintf(fp,"%f %f %f\n",vx,vy,vz);
              }
              else         fprintf(fp,"0.0 0.0 0.0\n");
              //fprintf(fp,"%f\n",(STLList[i]->Area[j]));
          }
      }
  }

  if(dump_what & DUMP_STRESSCOMPONENTS)
  {
     //write x y z component
     fprintf(fp,"VECTORS force float\n");
     for (int i = 0; i < STLList_len; i++)
     {
         double **f_tri = STLList[i]->f_tri;
         double * Area = STLList[i]->Area;
         for (int j=0;j<STLList[i]->nTri;j++)
            fprintf(fp, "%f\t%f\t%f\n",f_tri[j][0]/Area[j], f_tri[j][1]/Area[j],f_tri[j][2]/Area[j]);
      }
   }
  //footer not needed
  //if would be needed, would do like in dump stl
}
