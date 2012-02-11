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

// contributing author for rotation option: Evan Smuts (U Cape Town)

#ifndef LMP_STLTRI_H
#define LMP_STLTRI_H

#include "lammps.h"
#include "pointers.h"

namespace LAMMPS_NS {

class STLtri: protected Pointers {
 public:
  STLtri(LAMMPS *lmp);
  ~STLtri();

  void grow_arrays();
  void initMove(double);
  void initConveyor();
  void initRotation(double);
  void pointToVec();
  void vecToPoint();
  void getTriBoundBox(int, double *, double *,double);
  void getMeshBoundBox(double *xmin, double *xmax);
  void before_rebuild();
  bool are_coplanar_neighs(int i1,int i2);
  bool is_planar();
  void normal_vec(double *vec,int i);

  inline int mesh_moving(){return (movingMesh | conveyor | rotation);}
  inline double const* const* const* vnode() { return v_node; }

  int is_on_surface(double *pos);
  int is_in_tri(double *pos,int i);
  void generate_random(double *,class RanPark* random);

  int nTri,nTriMax;

  int movingMesh;  
  int conveyor;    
  double conv_vel[3];

  int rotation;
  double rot_axis[3];
  double rot_origin[3];
  double rot_omega;

  double Temp; //temperature of the mesh

  double skinSafetyFactor; 

  int* EDGE_INACTIVE;
  int* CORNER_INACTIVE;

  double Area_total;

   #define VECPERTRI 11

   //these are the triangles read from the STL file
   //index1:facet(triangle) index2:x y z

   double ***node;        
   double ***v_node;      
   double ***node_lastRe; 
   double **facenormal;
   double **f_tri;        
   double **fn_fshear;    
   double *Area;          
   double *wear;          
   double *wear_step;          

   double **cK;
   double ***ogK;
   double **ogKlen;
   
   double ***oKO; 
   double *rK;

   int **neighfaces; 
   int  *contactInactive;

  //the data is stored in the arrays above, the x array below have pointers to the data
  
  double **x;         
  double **xoriginal; 
  
  double **v;
  double **f;
  double *rmass;  
  int xvf_len;    
  int xvf_lenMax; 

 private:
  bool alreadyInitialized;

}; //end class

}

#endif

