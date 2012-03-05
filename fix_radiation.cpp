#include "fix_radiation.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "atom.h"
#include "domain.h"
#include "group.h"
#include "force.h"
#include "comm.h"
#include "update.h"
#include "error.h"
#include "modify.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "pair_gran.h"
#include "math_extra.h"
#include "fix_property_global.h"
#include "fix_property_atom.h"
#include "fix_scalar_transport_equation.h"
#include "mech_param_gran.h"
#include "respa.h"
#include "compute_pair_gran_local.h"

using namespace LAMMPS_NS;

FixRadiation::FixRadiation(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
    if (narg < 9) error->all("Illegal fix heat/gran command, not enough arguments");
    
    fix_temp = fix_heatFlux = fix_heatSource = NULL;
    fix_conductivity = NULL;
    
    conductivity = NULL;
    
    sx = atof(arg[3]);
    sy = atof(arg[4]);
    sz = atof(arg[5]);
    
    sr = atof(arg[6]);
    wavelength = atof(arg[7]);
    intensity = atof(arg[8]);
    
    ss = M_PI * pow(sr+sr,2);
    sp = ss*1*pow(TEMP,4)*BOLTS;
    
    pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
}

void FixRadiation::post_force(int a)
{  
  int *ilist,*jlist,*numneigh,**firstneigh;
  double x,y,z,radj,radi,nx,ny,nz,Cp,shapef;
  
  int i,j,ii,jj,inum,jnum;  
  
  int newton_pair = force->newton_pair;
  
  if (strcmp(force->pair_style,"hybrid")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid/overlay");
  
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
  
  updatePtrs();
  
  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
      
      i = ilist[ii];
      
      if (!(mask[i] & groupbit)) continue;
      
      x = pos[i][0]; 
      y = pos[i][1];
      z = pos[i][2];
      radi = radius[i];
      jlist = firstneigh[i];
      jnum = numneigh[i];
      
      //cout<<radi<<" "<<sr<<endl;
      
      Cp = conductivity[type[i]-1];
      //cout<<"Conductivity: "<<conductivity[type[i]-1]<<endl;
      
      double distc = dist(sx,sy,sz,x,y,z);
      //cout<<"Distance"<<distc<<endl;
      if(distc <= 10*radi)
      {
      shapef = procedeCalc(radi, sr, distc);
      //cout<<"Shape factor: "<<shapef<<endl;
      double T1 = Temp[i];
      //cout<<"Old temp: "<<T1<<endl;
      double m = rmass[i];
      //cout<<"Mass: "<<m<<endl;
      double T2 = T1 + (sp*shapef)/(m*Cp);
      
      heatFlux[i] = T2;
      //cout<<"Old temp: "<<T1<<", old TempFlux: "<<heatFlux[i]<<" new temp: "<<T2<<endl;
      }
  }      
}

double FixRadiation::dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double delx = x2-x1; double dely=y2-y1; double delz=z2-z1;
    return delx*delx + dely*dely + delz*delz;
}

FixRadiation::~FixRadiation()
{
    
}

void FixRadiation::updatePtrs()
{
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
}

void FixRadiation::init()
{
      if (!atom->radius_flag || !atom->rmass_flag) error->all("Please use a granular atom style for fix heat/gran");

      if(!force->pair_match("gran", 0)) error->all("Please use a granular pair style for fix heat/gran");
      
        pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
        
        fix_temp = static_cast<FixPropertyAtom*>(modify->find_fix_property("Temp","property/atom","scalar",0,0));
        fix_heatFlux = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatFlux","property/atom","scalar",0,0));
        fix_heatSource = static_cast<FixPropertyAtom*>(modify->find_fix_property("heatSource","property/atom","scalar",0,0));

        int max_type = pair_gran->mpg->max_type();

        if(conductivity) delete []conductivity;
        conductivity = new double[max_type];
        fix_conductivity = static_cast<FixPropertyGlobal*>(modify->find_fix_property("thermalConductivity","property/global","peratomtype",max_type,0));

        // pre-calculate conductivity for possible contact material combinations
        for(int i=1;i< max_type+1; i++)
            for(int j=1;j<max_type+1;j++)
            {
                conductivity[i-1] = fix_conductivity->compute_vector(i-1);
                if(conductivity[i-1] < 0.) error->all("Fix heat/gran: Thermal conductivity must not be < 0");
            }
      
      //cout<<"Calling updatePtrs()"<<endl;
      updatePtrs();
        
      cout<<"INIT finished"<<endl;
}


//Mask to provide to LAMMPS that we're not doing shit. Type is FORCE. Maybe change...
int FixRadiation::setmask()
{
        int mask = 0;
        mask |= POST_FORCE;
        return mask;
}

//To extend to a third D integral. (integral depending only on one dimension)
double third_D_coef(double d, double sind)
{
	return 2 * M_PI * d * sind;
}

//returns the squared norm of the vector vect (2 Dimensions)
double sqr_vect_norm(double * vect)
{
	return pow(vect[0], 2) + pow(vect[1], 2);
}

//Calculates the contribution of one element in the integral
double integrate(struct param * par)
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
	double scala = par->cosa * par->R[0] + par->sina * par->R[1];
	double scalb = (-par->cosb) * par->R[0] + (-par->sinb) * par->R[1];

	if(scala > 0)
		test = true;

	if(scalb > 0)
		test = true;

	return test;
}

double FixRadiation::procedeCalc(double radA, double radB, double dist)
{
	int ia = 0;
	int ib = 0;
	double integral = 0;

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
        return par.integral;
}