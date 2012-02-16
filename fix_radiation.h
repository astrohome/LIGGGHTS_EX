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

#include "fix.h"

namespace LAMMPS_NS {

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
};

}

#endif
#endif
