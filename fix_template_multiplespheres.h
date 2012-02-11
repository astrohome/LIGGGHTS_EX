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

FixStyle(particletemplate/multiplespheres,FixTemplateMultiplespheres)

#else

#ifndef LMP_FIX_TEMPLATE_MULTIPLESPHERES_H
#define LMP_FIX_TEMPLATE_MULTIPLESPHERES_H

#include "fix.h"
#include "fix_template_sphere.h"

namespace LAMMPS_NS {

class FixTemplateMultiplespheres : public FixTemplateSphere {
 public:
  FixTemplateMultiplespheres(class LAMMPS *, int, char **);
  ~FixTemplateMultiplespheres();
  virtual void post_create();
  double max_r_bound();
  double max_rad();
  int number_spheres();

  // single insertion
  virtual void randomize_single();

  // multi insertion
  void init_ptilist(int);
  void randomize_ptilist(int ,int );

  virtual void finalize_insertion() {}

 protected:

  // template calculations
  virtual void calc_bounding_sphere();
  virtual void calc_center_of_mass();

  // sqr distance from x_sphere[j] to xtest
  double dist_sqr(int j,double *xtest);

  // generate random point in bbox
  void generate_xtry(double *xtry);

  // number of spheres in template
  int nspheres;

  // coords of each sphere with respect to center of mass
  double **x_sphere;

  // radius of each sphere
  double *r_sphere;

  // scale factor if read from a file
  double scale_fact;

  // bounding box
  double x_min[3], x_max[3];

  // bounding sphere - radius and coordinates with respect to com
  double r_bound;
  double x_bound[3];

  // radius of sphere with equal volume
  double r_equiv;

  // number of tries for mc
  int ntry;
};

}

#endif
#endif

