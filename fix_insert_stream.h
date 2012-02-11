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

FixStyle(insert/stream,FixInsertStream)

#else

#ifndef LMP_FIX_INSERT_STREAM_H
#define LMP_FIX_INSERT_STREAM_H

#include "fix_insert.h"

namespace LAMMPS_NS {

class FixInsertStream : public FixInsert {
 public:

  FixInsertStream(class LAMMPS *, int, char **);
  ~FixInsertStream();
  void post_create();
  virtual int setmask();
  virtual void init();
  void init_defaults();
  virtual void end_of_step();

 protected:

  virtual void calc_insertion_properties();

  void pre_insert();

  int overlap(int);
  inline void generate_random(double *pos, double rad);

  void x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this);
  virtual void finalize_insertion(int);

  // additional insertion settings
  int duration;            //duration for insertion in time-steps

  // stuff for insertion region
  double normalvec[3];
  double extrude_length;
  double p_ref[3];         //reference point on face
  int face_style;
  double v_normal[3];      // insertion velocity projected on face

  // mesh face and bounding box of extruded face
  class FixMeshGran *ins_face;
  double ins_vol_xmin[3];
  double ins_vol_xmax[3];

  // non-mesh face
  double c_center[3], c_r; // non mesh face - currently only circle

  // per particle release step
  class FixPropertyAtom *fix_release;
  char *releasevar_name;

};

}

#endif
#endif
