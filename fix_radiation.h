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

#include "fix.h"
#include <iostream>
#include <math.h>

using namespace std;

namespace LAMMPS_NS {
    
struct param {
  float da;
  float db;
  float D;

  float cosa;
  float cosb;
  float sina;
  float sinb;

  float R [2];
  float integral;

};

class FixRadiation : public Fix {
 public:
  FixRadiation(class LAMMPS *, int, char **);
  ~FixRadiation();
  int setmask();
  void init(LAMMPS *);
  void post_force(int);
 private:
     int temp;
     NeighList *list;
     class PairGran *pair_gran;
     float dist(float,float,float,float,float,float);
     int procedeCalc(float, float, float);
     float sx, sy, sz, intensity, sr, wavelength;

};

}

#endif
#endif
