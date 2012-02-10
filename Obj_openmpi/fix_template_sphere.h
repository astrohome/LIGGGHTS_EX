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

#ifdef FIX_CLASS

FixStyle(particletemplate/sphere,FixTemplateSphere)

#else

#ifndef LMP_FIX_TEMPLATE_SPHERE_H
#define LMP_FIX_TEMPLATE_SPHERE_H

#include "fix.h"
#include "probability_distribution.h"

namespace LAMMPS_NS {

class FixTemplateSphere : public Fix {
 public:
  FixTemplateSphere(class LAMMPS *, int, char **);
  ~FixTemplateSphere();
  virtual void post_create(){}
  void write_restart(FILE *);
  void restart(char *);

  virtual int setmask();
  virtual double volexpect();           
  virtual double massexpect();          
  virtual double max_rad();
  virtual double max_r_bound();
  int number_spheres();
  int type();
  class Region *region();

  // ****single particle generation****
  virtual void randomize_single();              
  class ParticleToInsert *pti;

  // ****many particle generation****
  virtual void init_ptilist(int);
  virtual void delete_ptilist();
  virtual void randomize_ptilist(int,int);
  int n_pti_max;
  class ParticleToInsert **pti_list;

  virtual void finalize_insertion() {}

 protected:
  class Region *reg;
  class FixRegionVariable *reg_var;

  class RanPark *random;
  int seed;
  double PI;
  int iarg;

  // properties of particle template
  int nspheres;       
  double **x_sphere;   

  int atom_type;

  class LMP_PROBABILITY_NS::PDF *pdf_radius;   
  class LMP_PROBABILITY_NS::PDF *pdf_density;

  double volume_expect;
  double mass_expect;

};

}

#endif
#endif

