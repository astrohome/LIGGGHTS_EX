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
    int *ilist,*jlist,*numneigh,**firstneigh;
    double xtmp,ytmp,ztmp,radj,radi;
  
  int i,j,ii,jj,inum,jnum;
    
  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
    // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
      //Even don't ask me, what's going up here.
      i = ilist[ii];
      xtmp = x[i][0]; //Seems to be temperature. But why in 3D? :-/
      ytmp = x[i][1];
      ztmp = x[i][2];
      jlist = firstneigh[i]; //Neighbours???
      jnum = numneigh[i]; //No idea.
      radi = radius[i]; //Seems to be radius of current atom.
       for (jj = 0; jj < jnum; jj++) {
           j = jlist[jj];
           radj = radius[j];
           //printf("Radius?? A lot of them?? o_O %d\n",radj);
       }
      
      printf("Molecular MAN!! xtmp=%d, rad=%d\n",xtmp,radi);
  }
  
}

FixRadiation::~FixRadiation()
{
    
}

void FixRadiation::init(LAMMPS *lmp)
{

}


//Mask to provide to LAMMPS that we're not doing shit. Type is FORCE. Maybe change...
int FixRadiation::setmask()
{
        int mask = 0;
        mask |= POST_FORCE;
        return mask;
}
