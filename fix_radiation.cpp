#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_radiation.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"
#include "neighbor.h"
#include "atom.h"
#include "force.h"
#include "error.h"

using namespace LAMMPS_NS;

FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    if (narg < 4) error->all("Illegal fix heat/gran command, not enough arguments");
    
    printf("INITIALISATION!!!!!!!");
    pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
}

void FixRadiation::post_force(int a)
{    
    int *ilist,*numneigh,**firstneigh;
    double xtmp,ytmp,ztmp;
  
  int inum,ii,i;
    
    inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  atom->
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
    // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
      i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
      printf("Molecular MAN!! xtmp=%d\n",xtmp);
  }
  
}

FixRadiation::~FixRadiation()
{
    
}

void FixRadiation::init(LAMMPS *lmp)
{

}



int FixRadiation::setmask()
{
        int mask = 0;
        mask |= POST_FORCE;
        return mask;
}
