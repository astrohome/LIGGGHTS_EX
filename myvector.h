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

#ifndef LMP_MYVECTOR_H
#define LMP_MYVECTOR_H

#include<cmath>
#include "lammps.h"

namespace LAMMPS_NS {
//================================================
//SOME VERY SIMPLE VECTOR OPERATIONS
//================================================

inline double vectorMag3D(double *v)
{
  return (  std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2])  );
}

inline double vectorMag3DSquared(double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]  );
}

inline double vectorMag4D(double *v)
{
  return (  std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3])  );
}

inline double vectorMag4DSquared(double *v)
{
  return (  v[0]*v[0]+v[1]*v[1]+v[2]*v[2]+v[3]*v[3]  );
}

inline double vectorDot3D(double *v1,double *v2)
{
  return (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]);
}

inline void vectorCopy3D(double *from,double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
}

inline void vectorComponentMin3D(double *v1,double *v2,double *min)
{
    if(v1[0] > v2[0])
        min[0] = v2[0];
    else
        min[0] = v1[0];

    if(v1[1] > v2[1])
        min[1] = v2[1];
    else
        min[1] = v1[1];

    if(v1[2] > v2[2])
        min[2] = v2[2];
    else
        min[2] = v1[2];
}

inline void vectorComponentMax3D(double *v1,double *v2,double *max)
{
    if(v1[0] > v2[0])
        max[0] = v1[0];
    else
        max[0] = v2[0];

    if(v1[1] > v2[1])
        max[1] = v1[1];
    else
        max[1] = v2[1];

    if(v1[2] > v2[2])
        max[2] = v1[2];
    else
        max[2] = v2[2];
}

inline void vectorCopy4D(double *from,double *to)
{
  to[0]=from[0];
  to[1]=from[1];
  to[2]=from[2];
  to[3]=from[3];
}

inline void vectorScalarMult3D(double *v,double s)
{
  v[0]=s*v[0];
  v[1]=s*v[1];
  v[2]=s*v[2];
}

inline void vectorScalarMult3D(double *v,double s,double *result)
{
  result[0]=s*v[0];
  result[1]=s*v[1];
  result[2]=s*v[2];
}

inline void vectorScalarDiv3D(double *v,double s)
{
  v[0]=1./s*v[0];
  v[1]=1./s*v[1];
  v[2]=1./s*v[2];
}

inline void vectorNegate3D(double *v,double *result)
{
  result[0]=-v[0];
  result[1]=-v[1];
  result[2]=-v[2];
}

inline void vectorScalarDiv3D(double *v,double s,double *result)
{
  result[0]=1./s*v[0];
  result[1]=1./s*v[1];
  result[2]=1./s*v[2];
}

inline void vectorAdd3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]+v2[0];
  result[1]=v1[1]+v2[1];
  result[2]=v1[2]+v2[2];
}

inline void vectorSubtract3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[0]-v2[0];
  result[1]=v1[1]-v2[1];
  result[2]=v1[2]-v2[2];
}

inline void vectorCross3D(const double *v1,const double *v2, double *result)
{
  result[0]=v1[1]*v2[2]-v1[2]*v2[1];
  result[1]=v1[2]*v2[0]-v1[0]*v2[2];
  result[2]=v1[0]*v2[1]-v1[1]*v2[0];
}

inline void vectorZeroize3D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
}

inline void vectorZeroize4D(double *v)
{
  v[0]=0.;
  v[1]=0.;
  v[2]=0.;
  v[3]=0.;
}

inline void quatUnitize4D(double *q)
{
  q[0]=1.;
  q[1]=0.;
  q[2]=0.;
  q[3]=0.;
}

inline bool isUnitQuat4D(double *q)
{
    return
    (
        q[0] == 1. &&
        q[1] == 0. &&
        q[2] == 0. &&
        q[3] == 0.
    );
}

inline void normalize_bary(double *v)
{
  double mag = v[0]+v[1]+v[2];
  v[0]/=mag;
  v[1]/=mag;
  v[2]/=mag;
}

inline void vectorToBuf3D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
}

inline void bufToVector3D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
}

inline void vectorToBuf4D(double *vec,double *buf,int &m)
{
  buf[m++] = vec[0];
  buf[m++] = vec[1];
  buf[m++] = vec[2];
  buf[m++] = vec[3];
}

inline void bufToVector4D(double *vec,double *buf,int &m)
{
  vec[0] = buf[m++];
  vec[1] = buf[m++];
  vec[2] = buf[m++];
  vec[3] = buf[m++];
}

inline void printVec3D(FILE *out,char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e\n",name,vec[0],vec[1],vec[2]);
}

inline void printVec4D(FILE *out,char *name, double *vec)
{
    fprintf(out," vector %s: %e %e %e %e\n",name,vec[0],vec[1],vec[2],vec[3]);
}

inline void printMat33(FILE *out,char *name, double **mat)
{
    fprintf(out," matrix %s: %f %f %f\n",name,mat[0][0],mat[0][1],mat[0][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[1][0],mat[1][1],mat[1][2]);
    fprintf(out,"        %s: %f %f %f\n",name,mat[2][0],mat[2][1],mat[2][2]);
}

}

#endif
