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

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "memory.h"
#include "error.h"
#include "group.h"
#include "neighbor.h"
#include "fix_property_global.h"

using namespace LAMMPS_NS;

#define EPSILON 0.001
#define myAtof lmp->force->numeric 
#define myAtoi lmp->force->inumeric 

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::FixPropertyGlobal(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
    //Check args
    if (narg < 6) error->all("Illegal fix property/global command, not enough arguments");

    //Read args
    int n = strlen(arg[3]) + 1;
    variablename = new char[n];
    strcpy(variablename,arg[3]);

    if (strcmp(arg[4],"scalar") == 0) data_style = FIXPROPERTY_GLOBAL_SCALAR;
    else if (strcmp(arg[4],"vector") == 0) data_style = FIXPROPERTY_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"peratomtype") == 0 || strcmp(arg[4],"atomtype") == 0) data_style = FIXPROPERTY_GLOBAL_VECTOR;
    else if (strcmp(arg[4],"matrix") == 0) data_style = FIXPROPERTY_GLOBAL_MATRIX;
    else if (strcmp(arg[4],"peratomtypepair") == 0 || strcmp(arg[4],"atomtypepair") == 0) data_style = FIXPROPERTY_GLOBAL_MATRIX;
    else error->all("Unknown style for fix property/global. Valid styles are scalar or vector/peratomtype or matrix/peratomtypepair");

    int darg = 0;
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX) darg = 1;

    //assign values
    nvalues = narg - 5 - darg;
    
    values = (double*) memory->smalloc(nvalues*sizeof(double),"values");
    values_recomputed = (double*) memory->smalloc(nvalues*sizeof(double),"values");

    if(narg < 5+darg+nvalues) error->all("Illegal fix property/global command, not enough arguments");

    for (int j = 0; j < nvalues; j++)
        values[j] = myAtof(arg[5+darg+j]);

    if (data_style == FIXPROPERTY_GLOBAL_SCALAR)
        scalar_flag = 1;
    else if (data_style==FIXPROPERTY_GLOBAL_VECTOR) {
        vector_flag = 1;
        size_vector = nvalues;
    }
    else if (data_style == FIXPROPERTY_GLOBAL_MATRIX) {
        array_flag = 1;
        size_array_cols = myAtoi(arg[5]);
        if (fmod(static_cast<double>(nvalues),size_array_cols) != 0.)
          error->all("Error in fix property/global: The number of default values must thus be a multiple of the nCols.");
        size_array_rows = static_cast<int>(static_cast<double>(nvalues)/size_array_cols);
    }

    extvector=0; 

    //check if there is already a fix that tries to register a property with the same name
    for (int ifix = 0; ifix < modify->nfix; ifix++)
        if ((modify->fix[ifix]) && (strcmp(modify->fix[ifix]->style,style) == 0) && (strcmp(((FixPropertyGlobal*)(modify->fix[ifix]))->variablename,variablename)==0) )
            error->all("Error in fix property/global. There is already a fix that registers a variable of the same name");

    array = NULL;
    array_recomputed = NULL;
    if(data_style == FIXPROPERTY_GLOBAL_MATRIX)
    {
        array = (double**)memory->smalloc(size_array_rows*sizeof(double**),"FixPropGlob:array");
        array_recomputed = (double**)memory->smalloc(size_array_rows*sizeof(double**),"FixPropGlob:array_recomputed");
        for(int i = 0; i < size_array_rows; i++) array[i] = &values[i*size_array_cols];
        for(int i = 0; i < size_array_rows; i++) array_recomputed[i] = &values_recomputed[i*size_array_cols];
    }
}

/* ---------------------------------------------------------------------- */

FixPropertyGlobal::~FixPropertyGlobal()
{
  // delete locally stored arrays
  delete[] variablename;

  memory->sfree(values);
  memory->sfree(values_recomputed);

  if(array)            memory->sfree(array);
  if(array_recomputed) memory->sfree(array_recomputed);
}

/* ---------------------------------------------------------------------- */

Fix* FixPropertyGlobal::check_fix(const char *varname,const char *svmstyle,int len1,int len2,bool errflag)
{
    char errmsg[200];

    if(strcmp(varname,variablename) == 0)
    {
        if(strcmp(svmstyle,"scalar") == 0) len1 = 1;

        // check variable style
        if(
            (strcmp(svmstyle,"scalar") == 0 && data_style != FIXPROPERTY_GLOBAL_SCALAR) ||
            ((strcmp(svmstyle,"vector") == 0 || strcmp(svmstyle,"peratomtype") == 0) && data_style != FIXPROPERTY_GLOBAL_VECTOR) ||
            ((strcmp(svmstyle,"matrix") == 0 || strcmp(svmstyle,"peratomtypepair") == 0) && data_style != FIXPROPERTY_GLOBAL_MATRIX)
        )
        {
            if(errflag)
            {
                sprintf(errmsg,"%s",svmstyle);
                strcat(errmsg," style required for fix property/global variable ");
                strcat(errmsg,varname);
                error->all(errmsg);
            }
            else return NULL;
        }

        // check length
        if((nvalues < len1) && ((data_style != FIXPROPERTY_GLOBAL_MATRIX) || (data_style == FIXPROPERTY_GLOBAL_MATRIX) && (size_array_cols < len2)))
        {
            if(errflag)
            {
                sprintf(errmsg,"Length not sufficient for variable ");
                strcat(errmsg,varname);
                error->all(errmsg);
            }
            else return NULL;
        }

        // success
        return static_cast<Fix*>(this);
    }
    return NULL;
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::init()
{
    char errmsg[200];
    sprintf(errmsg,"Fix property/global: Length not sufficient for variable ");
    strcat(errmsg,variablename);
    if((strcmp(style,"property/atomtype") == 0 || strcmp(style,"property/peratomtype") == 0) && nvalues < atom->ntypes) error->all(errmsg);
    if((strcmp(style,"property/atomtypepair") == 0 || strcmp(style,"property/peratomtypepair") == 0) && nvalues < atom->ntypes*atom->ntypes) error->all(errmsg);
}

/* ---------------------------------------------------------------------- */

void FixPropertyGlobal::grow(int len1, int len2)
{
    if(data_style == FIXPROPERTY_GLOBAL_SCALAR) error->all("Can not grow global property of type scalar");
    else if(data_style == FIXPROPERTY_GLOBAL_VECTOR && len1 > nvalues)
    {
        values = (double*)memory->srealloc(values,len1*sizeof(double),"FixPropertyGlobal:values");
    }
    else if(data_style == FIXPROPERTY_GLOBAL_MATRIX && len1*len2 > nvalues)
    {
        values = (double*) memory->srealloc(values,len1*len2*sizeof(double),"FixPropertyGlobal:values");
        size_array_rows = len1;
        size_array_cols = len2;
        nvalues = len1*len2;
        array = (double**)memory->srealloc(array,size_array_rows*sizeof(double**),"FixPropGlob:array");
        for(int i = 0; i < size_array_rows; i++) array[i] = &values[i*size_array_cols];
    }
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_scalar()
{
  return values[0];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_vector(int i)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    return values[i];
}

void FixPropertyGlobal::vector_modify(int i,double val)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    values_recomputed[i] = val;
}

double FixPropertyGlobal::compute_vector_modified(int i)
{
    if (i>(nvalues-1))error->all("Trying to access vector in fix property/global, but index out of bounds");
    return values_recomputed[i];
}

/* ---------------------------------------------------------------------- */

double FixPropertyGlobal::compute_array(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    return array[i][j];
}

void FixPropertyGlobal::array_modify(int i, int j,double val) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    array_recomputed[i][j] = val;
}

double FixPropertyGlobal::compute_array_modified(int i, int j) //i is row, j is column
{
    if (i>(size_array_rows-1))error->all("Trying to access matrix in fix property/global, but row index out of bounds");
    if (j>(size_array_cols-1))error->all("Trying to access matrix in fix property/global, but column index out of bounds");

    return array_recomputed[i][j];
}

/* ---------------------------------------------------------------------- */

int FixPropertyGlobal::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double FixPropertyGlobal::memory_usage()
{
  double bytes = nvalues * sizeof(double);
  return bytes;
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

void FixPropertyGlobal::new_array(int l1,int l2)
{
    
    if (data_style == FIXPROPERTY_GLOBAL_MATRIX) error->all("Fix property/global: Can not allocate extra array for matrix style");
    array_flag = 1;
    size_array_rows = l1;
    size_array_cols = l2;
    nvalues = l1*l2;

    array = (double**)memory->create_2d_double_array(size_array_rows,size_array_cols,"FixPropGlob:array");
    array_recomputed = (double**)memory->create_2d_double_array(size_array_rows,size_array_cols,"FixPropGlob:array_recomputed");
}
