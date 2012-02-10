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

#ifdef FIX_CLASS

FixStyle(mesh/gran,FixMeshGran)

#else

#ifndef LMP_FIX_MESH_GRAN_H
#define LMP_FIX_MESH_GRAN_H

#include "fix.h"
#include "input.h"
#include "math.h"
#include "myvector.h"
#include "stl_tri.h"

namespace LAMMPS_NS {

class FixMeshGran : public Fix {
  friend class FixMoveTri;
  friend class FixTriNeighlist;
  friend class InputSTL;
  friend class InputMeshTri;
  friend class DumpSTL;
  friend class DumpMesh;
  friend class FixWallGran;

 public:
  FixMeshGran(class LAMMPS *, int, char **);
  ~FixMeshGran();
  virtual int setmask();
  void write_restart(FILE *);
  void restart(char *);
  virtual int write_restart_sub(FILE * fp,int n){return n;}
  virtual void restart_sub(char *){}

  virtual void final_integrate(){}
  virtual void init(){};
  virtual void post_integrate(){}
  int is_moving() {return STLdata->movingMesh;}
  inline double curvature(){return curvature_;}

  int atom_type_wall;//atom type that is assigned to the wall (to derive the mechanical properties) and thus to all pseudo wall particles
  double Temp_mesh; //wall temperature

  //these are the triangles read from the STL file
  //index1:facet(triangle) index2:x y z
  class STLtri *STLdata;
  double ***node;
  int nTri;

  double *p_ref;

 protected:
  double force_total[3];
  double torque_total[3];

  // assume surfaces as one curved face up to an angle of phi = acos(curvature)
  double curvature_;

  virtual void add_particle_contribution(double*,double*,int,double*,int,int){}

  bool analyseStress;
  int iarg;
  double scale_fact,*off_fact, *rot_angle; 

  virtual int n_children(){return 0;}
  virtual void children_write(FILE* fp) {}
  virtual void children_restart(double *){}

 private:
  
  int* EDGE_INACTIVE;
  int* CORNER_INACTIVE;

  void calcTriCharacteristics(int nTri,double ***node,double **cK,double ***ogK,double **ogKlen,double ***oKO,double *rK,double *Area,double &Area_total,double **facenormal,int **neighfaces,int *contactInactive);
  int get_max_index_sharedcorner(int iTri,int &nPrev,int *prevFaces,double *node2check,double ***node,double *rK,int **neighfaces);
}; //end class

}

#endif
#endif
