/* 
 * File:   fix_radiation.h
 * Author: felix
 *
 * Created on February 10, 2012, 4:50 PM
 */

#ifdef FIX_CLASS

FixStyle(radiation,FixRadiation)

#else

#ifndef FIX_RADIATION_H
#define	FIX_RADIATION_H

#include "fix.h"
namespace LAMMPS_NS {

class FixRadiation : public Fix {
 public:

  FixRadiation(class LAMMPS *, int, char **);
  ~FixRadiation();
  int setmask();
  void init();
  
 private:

};

}

#endif
#endif

