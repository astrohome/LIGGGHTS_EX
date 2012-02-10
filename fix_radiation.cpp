#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_radiation.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
      if (narg < 5) error->all("Illegal fix adapt command");
        nevery = atoi(arg[3]);
      if (nevery <= 0) error->all("Illegal fix adapt command");
}

FixRadiation::~FixRadiation()
{

}

FixRadiation::init()
{
    NeighList *list = neighbor->lists->listfull;
    for(int i=0;i<sizeof(list);i++)
    {
        
    }
    }
}