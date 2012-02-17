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
    double x,y,z,radj,radi;
  
  int i,j,ii,jj,inum,jnum;
    
  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **pos = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;
    // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
      //Even don't ask me, what's going up here.
      i = ilist[ii];
      x = pos[i][0]; //AAAAA, it's position!!! :-/
      y = pos[i][1];
      z = pos[i][2];
      jlist = firstneigh[i]; //Neighbours - this is NEIGHBOURS.
      jnum = numneigh[i]; //Count of neighbours.
      radi = radius[i]; //Seems to be radius of current atom.
       for (jj = 0; jj < jnum; jj++) {
           j = jlist[jj];
           radj = radius[j];
           //printf("Radius?? A lot of them?? o_O %d\n",radj);
       }
      //!!!!! ATTENTION!!!!! EPIC FAIL, float is %f in printf!!!
      printf("Molecular MAN!! (x,y,z)=(%f,%f,%f), rad=%f, neighbours count=%d\n",x,y,z,radi,jnum);
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
