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

#ifndef LMP_MY_MPI_H
#define LMP_MY_MPI_H

#include "mpi.h"

/* ---------------------------------------------------------------------- */

// a poor man's inline MPI wrappers for LIGGGHTS

/* ---------------------------------------------------------------------- */

namespace MyMPI {

  inline void My_MPI_Sum_Vector(double *vector,int,MPI_Comm comm);
  inline void My_MPI_Sum_Scalar(double &scalar,MPI_Comm comm);
  inline void My_MPI_Sum_Scalar(double &scalar,double &scalar_all,MPI_Comm comm);

  inline void My_MPI_Sum_Vector(int *vector,int,MPI_Comm comm);
  inline void My_MPI_Sum_Scalar(int &scalar,MPI_Comm comm);
  inline void My_MPI_Sum_Scalar(int &scalar,int &scalar_all,MPI_Comm comm);

  inline void My_MPI_Min_Scalar(double &scalar,MPI_Comm comm);
  inline void My_MPI_Min_Scalar(int &scalar,MPI_Comm comm);

  inline void My_MPI_Max_Scalar(double &scalar,MPI_Comm comm);
  inline void My_MPI_Max_Scalar(int &scalar,MPI_Comm comm);

};

inline void MyMPI::My_MPI_Sum_Vector(double *vector,int len, MPI_Comm comm)
{
    double *vector_all = new double [len];
    MPI_Allreduce(vector,vector_all,len,MPI_DOUBLE,MPI_SUM,comm);
    for(int i = 0; i < len; i++) vector[i] = vector_all[i];
    delete []vector_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Sum_Scalar(double &scalar,MPI_Comm comm)
{
    double scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_DOUBLE,MPI_SUM,comm);
    scalar = scalar_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Sum_Scalar(double &scalar,double &scalar_all,MPI_Comm comm)
{
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_DOUBLE,MPI_SUM,comm);
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Sum_Vector(int *vector,int len,MPI_Comm comm)
{
    int *vector_all = new int [len];
    MPI_Allreduce(vector,vector_all,len,MPI_INT,MPI_SUM,comm);
    for(int i = 0; i < len; i++) vector[i] = vector_all[i];
    delete []vector_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Sum_Scalar(int &scalar,MPI_Comm comm)
{
    int scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_INT,MPI_SUM,comm);
    scalar = scalar_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Sum_Scalar(int &scalar,int &scalar_all,MPI_Comm comm)
{
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_INT,MPI_SUM,comm);
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Min_Scalar(double &scalar,MPI_Comm comm)
{
    double scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_DOUBLE,MPI_MIN,comm);
    scalar = scalar_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Max_Scalar(double &scalar,MPI_Comm comm)
{
    double scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_DOUBLE,MPI_MAX,comm);
    scalar = scalar_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Min_Scalar(int &scalar,MPI_Comm comm)
{
    int scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_INT,MPI_MIN,comm);
    scalar = scalar_all;
}

/* ---------------------------------------------------------------------- */

inline void MyMPI::My_MPI_Max_Scalar(int &scalar,MPI_Comm comm)
{
    int scalar_all;
    MPI_Allreduce(&scalar,&scalar_all,1,MPI_INT,MPI_MAX,comm);
    scalar = scalar_all;
}

#endif
