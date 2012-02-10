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

#ifndef LMP_MATH_EXTRA_LIGGGHTS_H
#define LMP_MATH_EXTRA_LIGGGHTS_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "error.h"
#include "myvector.h"
#include "math_extra.h"

#define TOLERANCE_ORTHO 1e-10

namespace MathExtraLiggghts {

  inline void col_times3(const double m[3][3],const double *v, double *ans);

  inline double mdet(const double m[3][3],LAMMPS_NS::Error *error);

  //cubic root approx
  inline double cbrt_5d(double d);
  inline double cbrta_halleyd(const double a, const double R);
  inline double halley_cbrt1d(double d);

  inline int min(int a,int b);
  inline int max(int a,int b);
  inline int abs(int a);
  inline double min(double a,double b);
  inline double max(double a,double b);
  inline double abs(double a);

  inline void matrix_invert_4x4_special(double matrix[4][4]);
  inline int is_inside_tet(double *pos,double invmatrix[4][4]);

  inline void local_coosys_to_cartesian(double *global,double *local, double *ex_local, double *ey_local, double *ez_local);
  inline void cartesian_coosys_to_local(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error);
  inline void cartesian_coosys_to_local_orthogonal(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error);

  // quaternion oeprations
  inline void qconjugate(double *q, double *qc);
  inline void quat_from_vec(const double *v, double *q);
  inline void vec_from_quat(const double *q, double *v);
  inline void vec_quat_rotate(double *vec, double *quat, double *result);
  inline void vec_quat_rotate(double *vec, double *quat);
};

/* ----------------------------------------------------------------------
   matrix  times col vector 
------------------------------------------------------------------------- */

void MathExtraLiggghts::col_times3(const double m[3][3],const double *v, double *ans)
{
  ans[0] = m[0][0]*v[0]+v[1]*m[1][0]+v[2]*m[2][0];
  ans[1] = v[0]*m[0][1]+m[1][1]*v[1]+v[2]*m[2][1];
  ans[2] = v[0]*m[0][2]+v[1]*m[1][2]+m[2][2]*v[2];
}

/* ----------------------------------------------------------------------

Matrix determinant
------------------------------------------------------------------------- */

double MathExtraLiggghts::mdet(const double m[3][3],LAMMPS_NS::Error *error)
{
    return ( -m[0][2]*m[1][1]*m[2][0] + m[0][1]*m[1][2]*m[2][0] + m[0][2]*m[1][0]*m[2][1] - m[0][0]*m[1][2]*m[2][1] - m[0][1]*m[1][0]*m[2][2] + m[0][0]*m[1][1]*m[2][2] );

}

/* ----------------------------------------------------------------------
   Cubic root approx. 
------------------------------------------------------------------------- */

inline double MathExtraLiggghts::cbrt_5d(double d)
{
   const unsigned int B1 = 715094163;
   double t = 0.0;
   unsigned int* pt = (unsigned int*) &t;
   unsigned int* px = (unsigned int*) &d;
   pt[1]=px[1]/3+B1;
   return t;
}

inline double MathExtraLiggghts::cbrta_halleyd(const double a, const double R)
{
	const double a3 = a*a*a;
    const double b= a * (a3 + R + R) / (a3 + a3 + R);
	return b;
}

// cube root approximation using 1 iteration of Halley's method (double)
inline double MathExtraLiggghts::halley_cbrt1d(double d)
{
	double a = cbrt_5d(d);
	return cbrta_halleyd(a, d);
}

  int MathExtraLiggghts::min(int a,int b) { if (a<b) return a; return b;}
  int MathExtraLiggghts::max(int a,int b) { if (a>b) return a; return b;}
  int MathExtraLiggghts::abs(int a) { if (a>0) return a; return -a;}
  double MathExtraLiggghts::min(double a,double b) { if (a<b) return a; return b;}
  double MathExtraLiggghts::max(double a,double b) { if (a>b) return a; return b;}
  double MathExtraLiggghts::abs(double a) { if (a>0.) return a; return -a;}

/*----------------------------------------------------------------------
   inverts a special 4x4 matrix that looks like this

   m11 m12 m13 m14
   m21 m22 m23 m24
   m31 m32 m33 m34
   1   1   1   1
------------------------------------------------------------------------- */
inline void MathExtraLiggghts::matrix_invert_4x4_special(double matrix[4][4])
{
    double fac,invfac,v1x,v2x,v3x,v4x,v1y,v2y,v3y,v4y,v1z,v2z,v3z,v4z;
    v1x = matrix[0][0]; v1y = matrix[1][0]; v1z = matrix[2][0];
    v2x = matrix[0][1]; v2y = matrix[1][1]; v2z = matrix[2][1];
    v3x = matrix[0][2]; v3y = matrix[1][2]; v3z = matrix[2][2];
    v4x = matrix[0][3]; v4y = matrix[1][3]; v4z = matrix[2][3];

    fac = -v1x*v2z*v3y+v1x*v2y*v3z+v2z*v3y*v4x-v2y*v3z*v4x+v1x*v2z*v4y-
        v2z*v3x*v4y-v1x*v3z*v4y+v2x*v3z*v4y+v1z*
        (v2x*v3y-v3y*v4x+v2y*(-v3x+v4x)-v2x*v4y+v3x*v4y)-
        v1x*v2y*v4z+v2y*v3x*v4z+v1x*v3y*v4z-v2x*v3y*v4z+v1y*
        (v2z*v3x-v2x*v3z-v2z*v4x+v3z*v4x+v2x*v4z-v3x*v4z);

    invfac = 1./fac;

    matrix[0][0] = (-v3z*v4y+v2z*(-v3y+v4y)+v2y*(v3z-v4z)+v3y*v4z)*invfac;
    matrix[1][0] = (v1z*(v3y-v4y)+v3z*v4y-v3y*v4z+v1y*(-v3z+v4z))*invfac;
    matrix[2][0] = (-v2z*v4y+v1z*(-v2y+v4y)+v1y*(v2z-v4z)+v2y*v4z)*invfac;
    matrix[3][0] = (v1z*(v2y-v3y)+v2z*v3y-v2y*v3z+v1y*(-v2z+v3z))*invfac;

    matrix[0][1] = (v2z*(v3x-v4x)+v3z*v4x-v3x*v4z+v2x*(-v3z+v4z))*invfac;
    matrix[1][1] = (-v3z*v4x+v1z*(-v3x+v4x)+v1x*(v3z-v4z)+v3x*v4z)*invfac;
    matrix[2][1] = (v1z*(v2x-v4x)+v2z*v4x-v2x*v4z+v1x*(-v2z+v4z))*invfac;
    matrix[3][1] = (-v2z*v3x+v1z*(-v2x+v3x)+v1x*(v2z-v3z)+v2x*v3z)*invfac;

    matrix[0][2] = (-v3y*v4x+v2y*(-v3x+v4x)+v2x*(v3y-v4y)+v3x*v4y)*invfac;
    matrix[1][2] = (v1y*(v3x-v4x)+v3y*v4x-v3x*v4y+v1x*(-v3y+v4y))*invfac;
    matrix[2][2] = (-v2y*v4x+v1y*(-v2x+v4x)+v1x*(v2y-v4y)+v2x*v4y)*invfac;
    matrix[3][2] = (v1y*(v2x-v3x)+v2y*v3x-v2x*v3y+v1x*(-v2y+v3y))*invfac;

    matrix[0][2] = (v2z*v3y*v4x-v2y*v3z*v4x-v2z*v3x*v4y+v2x*v3z*v4y+v2y*v3x*v4z-v2x*v3y*v4z)*invfac;
    matrix[1][2] = (-v1z*v3y*v4x+v1y*v3z*v4x+v1z*v3x*v4y-v1x*v3z*v4y-v1y*v3x*v4z+v1x*v3y*v4z)*invfac;
    matrix[2][2] = (v1z*v2y*v4x-v1y*v2z*v4x-v1z*v2x*v4y+v1x*v2z*v4y+v1y*v2x*v4z-v1x*v2y*v4z)*invfac;
    matrix[3][2] = (-v1z*v2y*v3x+v1y*v2z*v3x+v1z*v2x*v3y-v1x*v2z*v3y-v1y*v2x*v3z+v1x*v2y*v3z)*invfac;
}

/*----------------------------------------------------------------------
   checks if a point is inside a tetrader, described by an inverse matrix
   This inverse matrix is the the inverse of a special 4x4 matrix that looks like this

   m11 m12 m13 m14
   m21 m22 m23 m24
   m31 m32 m33 m34
   1   1   1   1

   where m11,m21,m31 is vertex 1 etc.
------------------------------------------------------------------------- */

inline int MathExtraLiggghts::is_inside_tet(double *pos,double invmatrix[4][4])
{
    double result[4];

    result[0] = invmatrix[0][0] * pos[0] + invmatrix[0][1] * pos[1] + invmatrix[0][2] * pos[2] + invmatrix[0][3];
    result[1] = invmatrix[1][0] * pos[0] + invmatrix[1][1] * pos[1] + invmatrix[1][2] * pos[2] + invmatrix[1][3];
    result[2] = invmatrix[2][0] * pos[0] + invmatrix[2][1] * pos[1] + invmatrix[2][2] * pos[2] + invmatrix[2][3];
    result[3] = invmatrix[3][0] * pos[0] + invmatrix[3][1] * pos[1] + invmatrix[3][2] * pos[2] + invmatrix[3][3];

    if(max(result[0],max(result[1],max(result[2],result[3]))) > 1.0) return 0;
    return 1;
}

/*----------------------------------------------------------------------
   transform from local to global coords
------------------------------------------------------------------------- */

void MathExtraLiggghts::local_coosys_to_cartesian(double *global,double *local, double *ex_local, double *ey_local, double *ez_local)
{
    global[0] = local[0]*ex_local[0] + local[1]*ey_local[0] + local[2]*ez_local[0];
    global[1] = local[0]*ex_local[1] + local[1]*ey_local[1] + local[2]*ez_local[1];
    global[2] = local[0]*ex_local[2] + local[1]*ey_local[2] + local[2]*ez_local[2];
}

/*----------------------------------------------------------------------
   transform from global to local coords
------------------------------------------------------------------------- */

void MathExtraLiggghts::cartesian_coosys_to_local(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error)
{
  double M[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};
  double Mt[3][3] = {{0.,0.,0.},{0.,0.,0.},{0.,0.,0.}};

  // set up the matrix
  LAMMPS_NS::vectorCopy3D(ex_local,M[0]);
  LAMMPS_NS::vectorCopy3D(ey_local,M[1]);
  LAMMPS_NS::vectorCopy3D(ez_local,M[2]);
  MathExtra::transpose3(M,Mt);

  // solve
  MathExtra::mldivide3(Mt,global,local,error);
}

/*----------------------------------------------------------------------
   transform from global to local coords
   faster for orthogonal matrix
------------------------------------------------------------------------- */

void MathExtraLiggghts::cartesian_coosys_to_local_orthogonal(double *local,double *global, double *ex_local, double *ey_local, double *ez_local,LAMMPS_NS::Error *error)
{
  // check if orthogonal
  double dot1 = LAMMPS_NS::vectorDot3D(ex_local,ey_local);
  double dot2 = LAMMPS_NS::vectorDot3D(ey_local,ez_local);
  double dot3 = LAMMPS_NS::vectorDot3D(ez_local,ex_local);
  int flag = dot1 > TOLERANCE_ORTHO || dot2 > TOLERANCE_ORTHO || dot3 > TOLERANCE_ORTHO;
  if(flag) error->one("Insufficient accuracy: using MathExtraLiggghts::cartesian_coosys_to_local_orthogonal() for non-orthogonal coo-sys");

  // solve
  local[0] = global[0]*ex_local[0] + global[1]*ex_local[1] + global[2]*ex_local[2];
  local[1] = global[0]*ey_local[0] + global[1]*ey_local[1] + global[2]*ey_local[2];
  local[2] = global[0]*ez_local[0] + global[1]*ez_local[1] + global[2]*ez_local[2];
}

/* ----------------------------------------------------------------------
   conjugate of a quaternion: qc = conjugate of q
   assume q is of unit length
------------------------------------------------------------------------- */

void MathExtraLiggghts::qconjugate(double *q, double *qc)
{
  qc[0] = q[0];
  qc[1] = -q[1];
  qc[2] = -q[2];
  qc[3] = -q[3];
}

/* ----------------------------------------------------------------------
   construct quaternion4 from vector3
------------------------------------------------------------------------- */

void MathExtraLiggghts::quat_from_vec(const double *v, double *q)
{
  q[0] = 0.;
  q[1] = v[0];
  q[2] = v[1];
  q[3] = v[2];
}

/* ----------------------------------------------------------------------
   construct vector3 from quaternion4
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_from_quat(const double *q, double *v)
{
  v[0] = q[1];
  v[1] = q[2];
  v[2] = q[3];
}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_quat_rotate(double *vec, double *quat, double *result)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    MathExtra::multiply_quat_quat(quat,vecQ,temp);
    MathExtra::multiply_quat_quat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
}

/*----------------------------------------------------------------------
   rotoate vector by quaternion
------------------------------------------------------------------------- */

void MathExtraLiggghts::vec_quat_rotate(double *vec, double *quat)
{
    double vecQ[4], resultQ[4], quatC[4], temp[4], result[3];

    // construct quaternion (0,vec)
    quat_from_vec(vec,vecQ);

    // conjugate initial quaternion
    qconjugate(quat,quatC);

    // rotate by quaternion multiplications
    MathExtra::multiply_quat_quat(quat,vecQ,temp);
    MathExtra::multiply_quat_quat(temp,quatC,resultQ);

    // return result
    vec_from_quat(resultQ,result);
    LAMMPS_NS::vectorCopy3D(result,vec);
}

#endif

