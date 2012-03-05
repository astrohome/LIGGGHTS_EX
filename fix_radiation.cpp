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
#include "group.h"
#include "create_atoms.h"

using namespace LAMMPS_NS;

FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    if (narg < 9) error->all("Illegal fix heat/gran command, not enough arguments");
    
    sx = atof(arg[3]);
    sy = atof(arg[4]);
    sz = atof(arg[5]);
    
    sr = atof(arg[6]);
    wavelength = atof(arg[7]);
    intensity = atof(arg[8]);
    
    pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
}

void FixRadiation::post_force(int a)
{    
    int *ilist,*jlist,*numneigh,**firstneigh;
    double x,y,z,radj,radi,nx,ny,nz;
  
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
      
      i = ilist[ii];
      x = pos[i][0]; 
      y = pos[i][1];
      z = pos[i][2];
      radi = radius[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      cout<<radi<<" "<<sr<<endl;
      
      procedeCalc(radi, sr, dist(sx,sy,sz,x,y,z));
      
      //printf("Neighbours count=%d\n",jnum);
      
      /*for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;
      
          nx = pos[j][0];
          ny = pos[j][1];
          nz = pos[j][2];
          radj = radius[j];
          
      //printf("MOLECULAR MAN!! (x,y,z)=(%f,%f,%f), rad=%f\n",nx,ny,nz,radj);
                
      }*/
  }
  
}

float FixRadiation::dist(float x1, float y1, float z1, float x2, float y2, float z2)
{
    float delx = x2-x1; float dely=y2-y1; float delz=z2-z1;
    return delx*delx + dely*dely + delz*delz;
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

//To extend to a third D integral. (integral depending only on one dimension)
float third_D_coef(float d, float sind)
{
	return 2 * M_PI * d * sind;
}

//returns the squared norm of the vector vect (2 Dimensions)
float sqr_vect_norm(float * vect)
{
	return pow(vect[0], 2) + pow(vect[1], 2);
}

//Calculates the contribution of one element in the integral
float integrate(struct param * par)
{
	return (par->cosa * par->cosb / sqr_vect_norm(par->R));
}

void calculate_vector(struct param * par, int ia, int ib)
{
	//Calculates the value of the R vector
	par->R[0] = par->D - par->cosa - par->cosb;
	par->R[1] = - par->sina - par->sinb;
}

//This function tests if the two points can physically "see" each other,
//by testing the scalar product of the 2 vectors.
bool scalar_test(struct param * par)
{
	bool test = false;
	float scala = par->cosa * par->R[0] + par->sina * par->R[1];
	float scalb = (-par->cosb) * par->R[0] + (-par->sinb) * par->R[1];

	if(scala > 0)
		test = true;

	if(scalb > 0)
		test = true;

	return test;
}

int FixRadiation::procedeCalc(float radA, float radB, float dist)
{
	int ia = 0;
	int ib = 0;
	float integral = 0;

	//init parameter

	struct param par;

	//If the radius is not specified, take the default value
	if ( radA > 0 && radB > 0)
        {
		par.da = radA;
                par.db = radB;
        }
	else
		par.da = par.db = RADIUS;

	//If the distance is not specified, take the default value
	if( dist >= 0)
		par.D = dist;
	else
		par.D = DIST;

	par.integral = 0;

	int nb = NB_INT; //Number if integration steps

	//starting discrete integration
	for (ia=0;ia<nb;ia++)
	{
		integral = 0;

		// Calculates the vector da coordinates, according to the angle ia
		par.cosa = par.da * cos((M_PI * ia) / (2 * NB_INT));
		par.sina = par.da * sin((M_PI * ia) / (2 * NB_INT));

		for (ib=0;ib<nb;ib++)
		{

			// Calculates the vector db coordinates, according to the angle ib
			par.cosb = par.db * cos((M_PI * ib) / (2 * NB_INT));
			par.sinb = par.db * sin((M_PI * ib) / (2 * NB_INT));

			calculate_vector(&par, ia, ib);
			if(scalar_test(&par))
				integral = integral + integrate(&par) * third_D_coef(par.db, par.sinb);
		}
		integral = integral * third_D_coef(par.da, par.sina);
		par.integral = par.integral + integral;
	}

	par.integral = integral / (M_PI * 2 * M_PI * pow(par.da, 2)); //normalize integral
	cout 	<< "distance : " << par.D
			<< ", radius A : " << par.da << ", radius B : " << par.db
			<< ", shape factor = " << par.integral
			<< endl;

	return 0;
}
