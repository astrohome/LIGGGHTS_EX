#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_radiation.h"
#include "neighbor.h"
#include "neigh_list.h"

using namespace LAMMPS_NS;

void FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
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
    //Print for all atoms, its neighbors list

	NeighList **list = lmp->neighbor->lists;
    for(int i=0;i<list->inum;i++)
    {
    	printf("atom %d, neighbors :\n", list[i]->index);

    	for(int j=0;j<list[i]->inum;j++)
    	{
    		printf("	%d", list[i]->ilist[j]);
    	}
    }
}

