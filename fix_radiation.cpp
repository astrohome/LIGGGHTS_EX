#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_radiation.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"

using namespace LAMMPS_NS;

FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg < 5) error->all("Illegal fix adapt command");	
         nevery = atoi(arg[3]);
  if (nevery <= 0) error->all("Illegal fix adapt command");
         init(lmp);
}

FixRadiation::~FixRadiation()
{
    
}

void FixRadiation::init(LAMMPS *lmp)
{
    
}

int FixRadiation::setmask()
{
    return 0;
}
