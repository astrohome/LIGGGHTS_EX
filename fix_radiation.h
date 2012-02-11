/*
 * fix_radiation.h
 *
 *  Created on: Feb 11, 2012
 *      Author: jn
 */

#ifndef FIX_RADIATION_H_
#define FIX_RADIATION_H_

#include "fix.h"

namespace LAMMPS_NS {

class FixRadiation : public Fix {
public :
	FixRadiation(LAMMPS *, int , char **);
	void init(LAMMPS *);

};
}
#endif

#endif /* FIX_RADIATION_H_ */
