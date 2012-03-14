/* ----------------------------------------------------------------------
LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
Transfer Simulations

www.liggghts.com | www.cfdem.com
Christoph Kloss, christoph.kloss@cfdem.com

LIGGGHTS is based on LAMMPS
LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
http://lammps.sandia.gov, Sandia National Laboratories
Steve Plimpton, sjplimp@sandia.gov

Copyright (2003) Sandia Corporation. Under the terms of Contract
DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
certain rights in this software. This software is distributed under
the GNU General Public License.

See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_heat_gran.h"
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
#include <iostream>

using namespace LAMMPS_NS;
using namespace std;

#define SMALL 1e-8

/* ---------------------------------------------------------------------- */

FixHeatGran::FixHeatGran(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if ((!atom->radius_flag)||(!atom->rmass_flag)) error->all("Fix heat/gran needs per particle radius and mass");

  //check if a fix of this style already exists
  for (int i = 0; i < modify->nfix; i++)
      if (strcmp(modify->fix[i]->style,style) == 0) error->all("A fix of type heat/gran is already registered. Cannot have more than one");

  if (narg < 4) error->all("Illegal fix heat/gran command, not enough arguments");
  T0 = atof(arg[3]);
  
  nb_int = atof(arg[4]);

  fix_temp = fix_heatFlux = fix_heatSource = NULL;
  fix_conductivity = NULL;

  conductivity = NULL;

  peratom_flag = 1;              
  size_peratom_cols = 0;         
  peratom_freq = 1;
  time_depend = 1;

  scalar_flag = 1; 
  global_freq = 1; 

  cpl = NULL;

  area_correction = 0;
}

/* ---------------------------------------------------------------------- */

FixHeatGran::~FixHeatGran()
{
    
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_create()
{
    fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");

    if(!fix_ste)
    {
        char **newarg = new char*[15];
        newarg[0] = (char *) "ste_heattransfer";
        newarg[1] = group->names[igroup];
        newarg[2] = (char *) "transportequation/scalar";
        newarg[3] = (char *) "equation_id";
        newarg[4] = (char *) "heattransfer";
        newarg[5] = (char *) "quantity";
        newarg[6] = (char *) "Temp";
        newarg[7] = (char *) "default_value";
        newarg[8] = new char[30];
        sprintf(newarg[8],"%f",T0);
        newarg[9] = (char *) "flux_quantity";
        newarg[10] = (char *) "heatFlux";
        newarg[11] = (char *) "source_quantity";
        newarg[12] = (char *) "heatSource";
        newarg[13] = (char *) "capacity_quantity";
        newarg[14] = (char *) "thermalCapacity";
        modify->add_fix(15,newarg);
        delete [] newarg[8];
        delete [] newarg;
    }
}

/* ---------------------------------------------------------------------- */

int FixHeatGran::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}
/* ---------------------------------------------------------------------- */
void FixHeatGran::updatePtrs()
{
  Temp = fix_temp->vector_atom;
  vector_atom = Temp; 

  heatFlux = fix_heatFlux->vector_atom;
  heatSource = fix_heatSource->vector_atom;
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::init()
{
  const double *Y,*nu,*Y_orig;
  double expo,Yeff_ij,Yeff_orig_ij,ratio;
  Fix* ymo_fix;

  if (!atom->radius_flag || !atom->rmass_flag) error->all("Please use a granular atom style for fix heat/gran");

  if(!force->pair_match("gran", 0)) error->all("Please use a granular pair style for fix heat/gran");

  pair_gran = static_cast<PairGran*>(force->pair_match("gran", 0));
  history_flag = pair_gran->is_history();

  fix_ste = modify->find_fix_scalar_transport_equation("heattransfer");
  if(!fix_ste) error->all("Fix heat/gran needs a fix transportequation/scalar to work with");

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

  // calculate heat transfer correction

  ymo_fix = modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",0,0,false);
  if(ymo_fix)
  {
      area_correction = 1;
      if(screen && comm->me == 0) fprintf(screen,"Activating area correction for fix heat/gran\n");

      if(force->pair_match("gran/hooke",0)) expo = 1.;
      else if(force->pair_match("gran/hertz",0)) expo = 2./3.;
      else error->all("Fix heat/gran with area correction could not identify the granular pair style you are using, supported are hooke and hertz types");

      Y = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0))->get_values();
      nu = static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0))->get_values();
      Y_orig = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0))->get_values();

      // allocate a new array within youngsModulusOriginal
      static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0))->new_array(max_type,max_type);

      // feed deltan_ratio into this array
      for(int i = 1; i < max_type+1; i++)
      {
          for(int j = 1; j < max_type+1; j++)
          {
              Yeff_ij      = 1./((1.-pow(nu[i-1],2.))/Y[i-1]     +(1.-pow(nu[j-1],2.))/Y[j-1]);
              Yeff_orig_ij = 1./((1.-pow(nu[i-1],2.))/Y_orig[i-1]+(1.-pow(nu[j-1],2.))/Y_orig[j-1]);
              ratio = pow(Yeff_ij/Yeff_orig_ij,expo);
              
              static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0))->array_modify(i-1,j-1,ratio);
          }
      }

      // get reference to deltan_ratio
      deltan_ratio = static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulusOriginal","property/global","peratomtype",max_type,0))->get_array_modified();
  }

  updatePtrs();
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::post_force(int vflag)
{
    //template function for using touchflag or not
    if(history_flag == 0) post_force_eval<0>(vflag,0);
    if(history_flag == 1) post_force_eval<1>(vflag,0);
}

/* ---------------------------------------------------------------------- */

void FixHeatGran::cpl_evaluate(ComputePairGranLocal *caller)
{
    if(caller != cpl) error->all("Illegal situation in FixHeatGran::cpl_evaluate");
    if(history_flag == 0) post_force_eval<0>(0,1);
    if(history_flag == 1) post_force_eval<1>(0,1);
}

/* ---------------------------------------------------------------------- */

template <int HISTFLAG>
void FixHeatGran::post_force_eval(int vflag,int cpl_flag)
{
  double hc,contactArea,delta_n;
  int i,j,ii,jj,inum,jnum;
  double xtmp,ytmp,ztmp,delx,dely,delz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,tcoi,tcoj;
  
  double Q;  
  
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;

  int newton_pair = force->newton_pair;

  if (strcmp(force->pair_style,"hybrid")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid");
  if (strcmp(force->pair_style,"hybrid/overlay")==0) error->warning("Fix heat/gran implementation may not be valid for pair style hybrid/overlay");

  inum = pair_gran->list->inum;
  ilist = pair_gran->list->ilist;
  numneigh = pair_gran->list->numneigh;
  firstneigh = pair_gran->list->firstneigh;
  if(HISTFLAG) firsttouch = pair_gran->listgranhistory->firstneigh;

  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double **x = atom->x;
  int *type = atom->type;
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  updatePtrs();

  // loop over neighbors of my atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];
    
    Q = powerCalc(radi,Temp[i]);
    
    heatFlux[i] -= Q;
    int r;
    for(r=0;r<inum;r++)
    {
        j = ilist[r];
        if(i==j) continue;
        
        //cout<<"Test"<<endl;
        radj = radius[j];
        
        if(radi>radj)        
                heatFlux[j] += Q * procedeCalc(radj,radi,dist(xtmp, ytmp,ztmp,x[j][0],x[j][1],x[j][2]));
        else heatFlux[j] += Q * procedeCalc(radi,radj,dist(xtmp, ytmp,ztmp,x[j][0],x[j][1],x[j][2]));
    }
    
    
    if(HISTFLAG) touch = firsttouch[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      if (!(mask[i] & groupbit) && !(mask[j] & groupbit)) continue;

      if(!HISTFLAG)
      {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx*delx + dely*dely + delz*delz;
          radj = radius[j];
          radsum = radi + radj;
      }

      if (HISTFLAG && touch[jj] || !HISTFLAG && (rsq < radsum*radsum)) {  //contact
         
         if(HISTFLAG)
         {
              delx = xtmp - x[j][0];
              dely = ytmp - x[j][1];
              delz = ztmp - x[j][2];
              rsq = delx*delx + dely*dely + delz*delz;
              radj = radius[j];
              radsum = radi + radj;
         }

         r = sqrt(rsq);

         if(area_correction)
         {
             delta_n = radsum - r;
             delta_n *= deltan_ratio[type[i]-1][type[j]-1];
             r = radsum - delta_n;
         }

         contactArea = - M_PI/4 * ( (r-radi-radj)*(r+radi-radj)*(r-radi+radj)*(r+radi+radj) )/(r*r); //contact area of the two spheres

         tcoi = conductivity[type[i]-1];
         tcoj = conductivity[type[j]-1];
         if (tcoi < SMALL || tcoj < SMALL) hc = 0.;
         else hc = 4.*tcoi*tcoj/(tcoi+tcoj)*sqrt(contactArea);

         if(!cpl_flag)
         {
             heatFlux[i] += (Temp[j]-Temp[i])*hc;
             if (newton_pair || j < nlocal) heatFlux[j] += (Temp[i]-Temp[j])*hc;
         }

         if(cpl_flag && cpl) cpl->add_heat(i,j,(Temp[j]-Temp[i])*hc);
      }
    }
  }

  if(newton_pair) fix_heatFlux->do_reverse_comm();
}

/* ---------------------------------------------------------------------- */
double FixHeatGran::compute_scalar()
{
    return fix_ste->compute_scalar();
}

/* ----------------------------------------------------------------------
   register and unregister callback to compute
------------------------------------------------------------------------- */

void FixHeatGran::register_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != NULL) error->all("Fix heat/gran allows only one compute of type pair/local");
   cpl = ptr;
}

void FixHeatGran::unregister_compute_pair_local(ComputePairGranLocal *ptr)
{
   
   if(cpl != ptr) error->all("Illegal situation in PairGran::unregister_compute_pair_local");
   cpl = NULL;
}




double FixHeatGran::powerCalc(double radius, double temperature)
{
    double surf = M_PI * pow(2*radius,2);
    return surf*1*pow(temperature,4) * BOLTS;
}


double FixHeatGran::dist(double x1, double y1, double z1, double x2, double y2, double z2)
{
    double delx = x2-x1; double dely=y2-y1; double delz=z2-z1;
    return sqrt(delx*delx + dely*dely + delz*delz);
}


//To extend to a third D integral. (integral depending only on one dimension)
double FixHeatGran::third_D_coef(double d, double sind, int nb)
{
	return 2 * pow((M_PI * d), 2) * sind / ((double)nb);
}

//returns the squared norm of the vector vect (2 Dimensions)
double FixHeatGran::sqr_vect_norm(double vecta, double vectb)
{
	return pow(vecta, 2) + pow(vectb, 2);
}

//Calculates the contribution of one element in the integral
double FixHeatGran::integrate(struct param * par)
{
	double cosg1, cosg2;
	double normR;
	normR = sqrt(sqr_vect_norm(par->R[0], par->R[1]));
	cosg1 = (par->xda * par->R[0] + par->yda * par->R[1]) / (par->da * normR);
 	cosg2 = (par->xdb * par->R[0] + par->ydb * par->R[1]) / (par->db * normR);

 	if (normR == 0)
 		return 0;
 	else
 		return cosg1 * cosg2 / pow(normR, 2);
}

void FixHeatGran::calculate_vector(struct param * par, int ia, int ib)
{
	//Calculates the value of the R vector
	par->R[0] = par->D - par->xda - par->xdb;
	par->R[1] = - par->yda - par->ydb;

}

//This function tests if the two points can physically "see" each other,
//by testing the scalar product of the 2 vectors.
bool FixHeatGran::scalar_test(struct param * par)
{
	double scala = par->xda * par->R[0] + par->yda * par->R[1];
	double scalb = par->xdb * par->R[0] + par->ydb * par->R[1];

	if(scala > 0 && scalb > 0)
		return true;

	else
		return false;
}

double FixHeatGran::procedeCalc(double radA, double radB, double dist)
{
	int ia = 0;
	int ib = 0;
	double integral = 0;

	//init parameter

	struct param par;

	//If the radius is not specified, take the default value
	par.da = radA;
	par.db = radB;
	par.D = dist;

	par.integral = 0;

	//starting discrete integration
	for (ia = 0 ; ia < nb_int ; ia++)
	{
		integral = 0;

		// Calculates the vector da coordinates, according to the angle ia
		par.xda = par.da * cos((M_PI * ia) / (2 * nb_int));
		par.yda = par.da * sin((M_PI * ia) / (2 * nb_int));

		for (ib=0;ib<2*nb_int;ib++)
		{

			// Calculates the vector db coordinates, according to the angle ib
			par.xdb = par.db * cos(-M_PI / 2 + (M_PI * ib) / (2 * nb_int));
			par.ydb = par.db * sin((M_PI * ib) / (2 * nb_int));

			calculate_vector(&par, ia, ib);
			if(scalar_test(&par))
				integral = integral + integrate(&par) * third_D_coef(par.db, par.ydb, nb_int) / 2;
		}
		integral = integral * third_D_coef(par.da, par.yda, nb_int);
		par.integral += integral;
	}
        par.integral /= 4 * pow((M_PI * par.da ),2); 
	return par.integral*80000;
}
