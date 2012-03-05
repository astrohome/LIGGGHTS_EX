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

#define NB_INT 	100
#define RADIUS 	5
#define DIST 	10 * RADIUS
#define TEMP    600
#define BOLTS   1.3806488*pow(10,-23)

#include "fix.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace LAMMPS_NS {
    
struct param {
  double da;
  double db;
  double D;

  double cosa;
  double cosb;
  double sina;
  double sinb;

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
     int temp;
     NeighList *list;
     class PairGran *pair_gran;
     double dist(double,double,double,double,double,double);
     double procedeCalc(double, double, double);
     double sx, sy, sz, intensity, sr, wavelength, ss, sp;
     
     class FixPropertyGlobal* fix_conductivity;
     class FixPropertyAtom* fix_temp;
     class FixPropertyAtom* fix_heatSource;
     
     double *conductivity;
     double *heatSource; 
     double *heatFlux; 
     double *Temp; 
     class FixPropertyAtom* fix_heatFlux;

};

}

#endif
#endif
