/* 
 * File:   fix_radiation.h
 * Author: felix
 *
 * Created on February 14, 2012, 10:49 PM
 */

#ifdef FIX_CLASS

FixStyle(heat/radiation,FixRadiation)

#else

#ifndef LMP_FIX_RADIATION_H
#define	LMP_FIX_RADIATION_H

#define BOLTS   5.670400*pow((double)10,-8)

#include "fix.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace LAMMPS_NS {
    
struct param {
  double da;
  double db;
  double D;

  double xda;
  double xdb;
  double yda;
  double ydb;

  double R [2];
  double integral;

};

class FixRadiation : public Fix {
 public:
  
  FixRadiation(class LAMMPS *, int, char **);
  ~FixRadiation();
  int setmask();
  void init();
  void post_force(int);
  void updatePtrs();
  
 private:
   
     NeighList *list;
     class PairGran *pair_gran;
     double dist(double,double,double,double,double,double);
     double procedeCalc(double, double, double);
     double powerCalc(double, double);
     double nb_int;     //Integration precision
     
     class FixPropertyAtom* fix_temp;
     class FixPropertyAtom* fix_heatFlux;
     class FixPropertyAtom* fix_heatSource;
     class FixPropertyGlobal* fix_capacity;
     class FixScalarTransportEquation *fix_ste;
     
     double *capacity;
     double *heatSource; 
     double *heatFlux; 
     double *Temp;

};

}

#endif
#endif
