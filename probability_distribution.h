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

#ifndef LMP_PROBABILITY_DISTRIBUTION_H
#define LMP_PROBABILITY_DISTRIBUTION_H

#include "math.h"
#include "stdio.h"
#include "string.h"
#include "random_park.h"
#include "error.h"
#include "pointers.h"

enum{RANDOM_CONSTANT,RANDOM_UNIFORM,RANDOM_GAUSSIAN,RANDOM_LOGNORMAL};

namespace LMP_PROBABILITY_NS {

   #define MATH_E 2.718281828459

  class PDF
  {
    public:
      PDF(LAMMPS_NS::Error *error)
      {
          mu_ = sigma_ = min_ = max_ = 0.;
          mass_shift_ = 0;
          this->error = error;
      }
      ~PDF(){}

      inline int rand_style()
      {
          return rand_style_;
      }

      inline void set_min_max(double min,double max)
      {
          min_ = min;
          max_ = max;
      }

      inline void activate_mass_shift()
      {
          mass_shift_ = 1;
      }

      template<int RAND_STYLE> void set_params(double)
      {
          error->all("Faulty usage of Probability::set_params");
      }

      template<int RAND_STYLE> void set_params(double,double)
      {
          error->all("Faulty usage of Probability::set_params");
      }

      int rand_style_;

      double mu_,sigma_;
      double min_,max_;

      // if 1, pdf is shifted from number to mass based pdf
      int mass_shift_;

      LAMMPS_NS::Error *error;
  };

  inline double pdf_max(PDF *pdf)
  {
      return pdf->max_;
  }

  inline double pdf_min(PDF *pdf)
  {
      return pdf->min_;
  }

  template <int RAND_STYLE> inline double expectancy_value(PDF *pdf)
  {
      pdf->error->all("Faulty usage of Probability::expectancy");
      return 0.;
  }

  template <int RAND_STYLE> inline double rand_value(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
      pdf->error->all("Faulty usage of Probability::rand");
      return 0.;
  }

  //------------------------------------------------------------------------------
  // CONSTANT
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_CONSTANT>(double val)
  {
	 rand_style_ = RANDOM_CONSTANT;
	 mu_ = val;
	 set_min_max(mu_,mu_);
  }

  template<> inline double expectancy_value<RANDOM_CONSTANT>(PDF *pdf)
  {

     return pdf->mu_;
  }

  template<> inline double rand_value<RANDOM_CONSTANT>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     return pdf->mu_;
  }

  //------------------------------------------------------------------------------
  // UNIFORM
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_UNIFORM>(double min, double max)
  {
	 rand_style_ = RANDOM_UNIFORM;
	 set_min_max(min,max);
  }

  template<> inline double expectancy_value<RANDOM_UNIFORM>(PDF *pdf)
  {
     if(!pdf->mass_shift_) return 0.5 * (pdf->min_ + pdf->max_);
     else return ( pow((pdf->max_*pdf->max_*pdf->max_*pdf->max_- pdf->min_*pdf->min_*pdf->min_*pdf->min_)/(4.*(pdf->max_ - pdf->min_)),0.333333) );
  }

  template<> inline double rand_value<RANDOM_UNIFORM>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     double rn =  rp->uniform();
     if(pdf->mass_shift_) rn = rn*rn*rn;
     return (pdf->min_) + rn * (pdf->max_ - pdf->min_);
  }

  //------------------------------------------------------------------------------
  // GAUSSIAN
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_GAUSSIAN>(double mu, double sigma)
  {
	 rand_style_ = RANDOM_GAUSSIAN;
	 mu_ = mu;
	 sigma_ = sigma;

	 //set min-max to +- 3 sigma (99.73% of all values)
	 set_min_max(mu_-3.*sigma_, mu_+3.*sigma_);
  }

  template<> inline double expectancy_value<RANDOM_GAUSSIAN>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all("mass_shift_ not implemented for gaussian");
     return pdf->mu_;
  }

  template<> inline double rand_value<RANDOM_GAUSSIAN>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     if(pdf->mass_shift_) pdf->error->all("mass_shift_ not implemented for gaussian");
     double value;
     do
     {
         value = pdf->mu_ + rp->gaussian() * pdf->sigma_;
     } while (value < pdf->min_ || value > pdf->max_);
     return value;
  }

  //------------------------------------------------------------------------------
  // LOGNORMAL
  //------------------------------------------------------------------------------

  template<> inline void PDF::set_params<RANDOM_LOGNORMAL>(double mu, double sigma)
  {
      rand_style_ = RANDOM_LOGNORMAL;
	  mu_ = mu;
	  sigma_ = sigma;

	  //also here, take +- 3 sigma as min/max
	  //change in expectancy considered negligable
	  double min =  pow(MATH_E,mu_ - 3. * sigma_);
	  double max =  pow(MATH_E,mu_ + 3. * sigma_);
	  set_min_max(min, max);
  }

  template<> inline double expectancy_value<RANDOM_LOGNORMAL>(PDF *pdf)
  {
     if(pdf->mass_shift_) pdf->error->all("mass_shift_ not implemented for lognormal");
     return pow(MATH_E,pdf->mu_ + 0.5 * pdf->sigma_ * pdf->sigma_);
  }

  template<> inline double rand_value<RANDOM_LOGNORMAL>(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
     if(pdf->mass_shift_) pdf->error->all("mass_shift_ not implemented for lognormal");
     double value;
     do
     {
         value = pow(MATH_E,pdf->mu_ + rp->gaussian() * pdf->sigma_);
     } while (value < pdf->min_ || value > pdf->max_);
     return value;
  }

  //------------------------------------------------------------------------------
  // MASTER FUNCTIONS
  //------------------------------------------------------------------------------

  inline double expectancy(PDF *pdf)
  {
      if(pdf->rand_style_ == RANDOM_CONSTANT) return expectancy_value<RANDOM_CONSTANT>(pdf);
      else if(pdf->rand_style_ == RANDOM_UNIFORM) return expectancy_value<RANDOM_UNIFORM>(pdf);
      else if(pdf->rand_style_ == RANDOM_GAUSSIAN) return expectancy_value<RANDOM_GAUSSIAN>(pdf);
      else if(pdf->rand_style_ == RANDOM_LOGNORMAL) return expectancy_value<RANDOM_LOGNORMAL>(pdf);
      else pdf->error->all("Faulty implemantation in Probability::expectancy");
      return 0.;
  }

  inline double rand(PDF *pdf,LAMMPS_NS::RanPark *rp)
  {
      if(pdf->rand_style_ == RANDOM_CONSTANT) return rand_value<RANDOM_CONSTANT>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_UNIFORM) return rand_value<RANDOM_UNIFORM>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_GAUSSIAN) return rand_value<RANDOM_GAUSSIAN>(pdf,rp);
      else if(pdf->rand_style_ == RANDOM_LOGNORMAL) return rand_value<RANDOM_LOGNORMAL>(pdf,rp);
      else pdf->error->all("Faulty implemantation in Probability::rand");
      return 0.;
  }

};

#endif
