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

/* ----------------------------------------------------------------------
   Contributing authors for original version: Leo Silbert (SNL), Gary Grest (SNL)
------------------------------------------------------------------------- */
#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "fix_wall_gran_hooke.h"
#include "pair_gran_hooke.h"
#include "atom.h"
#include "domain.h"
#include "update.h"
#include "force.h"
#include "pair.h"
#include "modify.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "compute_pair_gran_local.h"

using namespace LAMMPS_NS;

#define BIG 1.0e20

#define MIN(A,B) ((A) < (B)) ? (A) : (B)
#define MAX(A,B) ((A) > (B)) ? (A) : (B)

/* ---------------------------------------------------------------------- */

FixWallGranHooke::FixWallGranHooke(LAMMPS *lmp, int narg, char **arg) :
  FixWallGranHookeHistory(lmp, narg, arg)
{

}

/* ---------------------------------------------------------------------- */

void FixWallGranHooke::compute_force(int ip, double deltan, double rsq,double meff_wall, double dx, double dy, double dz,double *vwall,double *c_history,double area_ratio)
{
  double r,vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wrmag;
  double wr1,wr2,wr3,damp,ccel,vtr1,vtr2,vtr3,vrel;
  double fn,fs,ft,fs1,fs2,fs3,fx,fy,fz,tor1,tor2,tor3,rinv,rsqinv;
  double kn, kt, gamman, gammat, xmu, rmu;

  double *f = atom->f[ip];
  double *torque = atom->torque[ip];
  double *v = atom->v[ip];
  double *omega = atom->omega[ip];
  double radius = atom->radius[ip];
  double mass = atom->rmass[ip];
  double cr = radius - 0.5*deltan;

  //get the parameters needed to resolve the contact
  //deltan > 0 needed in this function
  deriveContactModelParams(ip,-deltan,meff_wall,kn,kt,gamman,gammat,xmu,rmu);

  r = sqrt(rsq);
  rinv = 1.0/r;
  rsqinv = 1.0/rsq;

  // relative translational velocity

  vr1 = v[0] - vwall[0];
  vr2 = v[1] - vwall[1];
  vr3 = v[2] - vwall[2];

  // normal component

  vnnr = vr1*dx + vr2*dy + vr3*dz;
  vn1 = dx*vnnr * rsqinv;
  vn2 = dy*vnnr * rsqinv;
  vn3 = dz*vnnr * rsqinv;

  // tangential component

  vt1 = vr1 - vn1;
  vt2 = vr2 - vn2;
  vt3 = vr3 - vn3;

  // relative rotational velocity

  wr1 = cr*omega[0] * rinv;
  wr2 = cr*omega[1] * rinv;
  wr3 = cr*omega[2] * rinv;

  // normal forces = Hookian contact + normal velocity damping

  damp = gamman*vnnr*rsqinv;   
  ccel = kn*(radius-r)*rinv - damp;
  
  // relative velocities

  vtr1 = vt1 - (dz*wr2-dy*wr3);
  vtr2 = vt2 - (dx*wr3-dz*wr1);
  vtr3 = vt3 - (dy*wr1-dx*wr2);
  vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
  vrel = sqrt(vrel);

  // force normalization

  fn = xmu * fabs(ccel*r);
  fs = gammat*vrel;         
  if (vrel != 0.0) ft = MIN(fn,fs) / vrel;
  else ft = 0.0;

  // tangential force due to tangential velocity damping

  fs1 = -ft*vtr1;
  fs2 = -ft*vtr2;
  fs3 = -ft*vtr3;

  // forces & torques

  fx = dx*ccel + fs1;
  fy = dy*ccel + fs2;
  fz = dz*ccel + fs3;

  if(addflag)
  {
      f[0] += fx*area_ratio;
      f[1] += fy*area_ratio;
      f[2] += fz*area_ratio;
  }

  tor1 = rinv * (dy*fs3 - dz*fs2);
  tor2 = rinv * (dz*fs1 - dx*fs3);
  tor3 = rinv * (dx*fs2 - dy*fs1);

  if(rollingflag)
  {
	    wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
	    if (wrmag > 0.)
	    {
	        tor1 += rmu*kn*(radius-r)*wr1/wrmag;
            tor2 += rmu*kn*(radius-r)*wr2/wrmag;
            tor3 += rmu*kn*(radius-r)*wr3/wrmag;
	    }
  }

  if(addflag)
  {
      torque[0] -= cr*tor1*area_ratio;
      torque[1] -= cr*tor2*area_ratio;
      torque[2] -= cr*tor3*area_ratio;
  }
  else if(cwl) cwl->add_wall_2(ip,fx,fy,fz,tor1*area_ratio,tor2*area_ratio,tor3*area_ratio,c_history,rsq);
}

