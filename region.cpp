/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "math.h"
#include "stdlib.h"
#include "string.h"
#include "region.h"
#include "update.h"
#include "domain.h"
#include "lattice.h"
#include "error.h"
#include "random_park.h"
#include "myvector.h"
#include "mympi.h"
#include "comm.h"

using namespace LAMMPS_NS;

enum{NONE,VELOCITY,WIGGLE,ROTATE,VARIABLE};

/* ---------------------------------------------------------------------- */

Region::Region(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  n = strlen(arg[1]) + 1;
  style = new char[n];
  strcpy(style,arg[1]);

  time_origin = update->ntimestep;
}

/* ---------------------------------------------------------------------- */

Region::~Region()
{
  delete [] id;
  delete [] style;

  if (random) delete random;
}

/* ---------------------------------------------------------------------- */

void Region::init()
{
  dt = update->dt;
}

/* ----------------------------------------------------------------------
   parse optional parameters at end of region input line
------------------------------------------------------------------------- */

void Region::options(int narg, char **arg)
{
  if (narg < 0) error->all("Illegal region command");

  // option defaults

  interior = 1;
  scaleflag = 1;
  dynamic = NONE;

  seed = 3012211;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"units") == 0) {
      if (iarg+2 > narg) error->all("Illegal region command");
      if (strcmp(arg[iarg+1],"box") == 0) scaleflag = 0;
      else if (strcmp(arg[iarg+1],"lattice") == 0) scaleflag = 1;
      else error->all("Illegal region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"side") == 0) {
      if (iarg+2 > narg) error->all("Illegal region command");
      if (strcmp(arg[iarg+1],"in") == 0) interior = 1;
      else if (strcmp(arg[iarg+1],"out") == 0) interior = 0;
      else error->all("Illegal region command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"vel") == 0) {
      if (iarg+4 > narg) error->all("Illegal region command");
      vx = atof(arg[iarg+1]);
      vy = atof(arg[iarg+2]);
      vz = atof(arg[iarg+3]);
      dynamic = VELOCITY;
      iarg += 4;
    } else if (strcmp(arg[iarg],"wiggle") == 0) {
      if (iarg+5 > narg) error->all("Illegal region command");
      ax = atof(arg[iarg+1]);
      ay = atof(arg[iarg+2]);
      az = atof(arg[iarg+3]);
      period = atof(arg[iarg+4]);
      dynamic = WIGGLE;
      iarg += 5;
    } else if (strcmp(arg[iarg],"rotate") == 0) {
      if (iarg+8 > narg) error->all("Illegal region command");
      point[0] = atof(arg[iarg+1]);
      point[1] = atof(arg[iarg+2]);
      point[2] = atof(arg[iarg+3]);
      axis[0] = atof(arg[iarg+4]);
      axis[1] = atof(arg[iarg+5]);
      axis[2] = atof(arg[iarg+6]);
      period = atof(arg[iarg+7]);
      dynamic = ROTATE;
      iarg += 8;
    
    } else if (strcmp(arg[iarg],"seed") == 0) {
      if (iarg+2 > narg) error->all("Illegal region command");
      seed = atoi(arg[iarg+1]);
    }
    else error->all("Illegal region command");
  }

  random = new RanPark(lmp,seed);

  // error check

  if (dynamic &&
      (strcmp(style,"union") == 0 || strcmp(style,"intersect") == 0))
    error->all("Region union or intersect cannot be dynamic");

  // setup scaling

  if (scaleflag && domain->lattice == NULL)
    error->all("Use of region with undefined lattice");

  if (scaleflag) {
    xscale = domain->lattice->xlattice;
    yscale = domain->lattice->ylattice;
    zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;

  if (dynamic == VELOCITY) {
    vx *= xscale;
    vy *= yscale;
    vz *= zscale;
  } else if (dynamic == WIGGLE) {
    ax *= xscale;
    ay *= yscale;
    az *= zscale;
  } else if (dynamic == ROTATE) {
    point[0] *= xscale;
    point[1] *= yscale;
    point[2] *= zscale;
  }

  if (dynamic == WIGGLE || dynamic == ROTATE) {
    double PI = 4.0 * atan(1.0);
    omega_rotate = 2.0*PI / period;
  }

  // runit = unit vector along rotation axis

  if (dynamic == ROTATE) {
    double len = sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
    if (len == 0.0)
      error->all("Region cannot have 0 length rotation vector");
    runit[0] = axis[0]/len;
    runit[1] = axis[1]/len;
    runit[2] = axis[2]/len;
  }
}

/* ----------------------------------------------------------------------
   return 1 if region is dynamic, 0 if static
   only primitive regions define it here
   union/intersect regions have their own dynamic_check()
------------------------------------------------------------------------- */

int Region::dynamic_check()
{
  return dynamic;
}

/* ----------------------------------------------------------------------
   determine if point x,y,z is a match to region volume
   XOR computes 0 if 2 args are the same, 1 if different
   note that inside() returns 1 for points on surface of region
   thus point on surface of exterior region will not match
   if region is dynamic, apply inverse of region change to x
------------------------------------------------------------------------- */

int Region::match(double x, double y, double z)
{
  double a[3],b[3],c[3],d[3];

  if (dynamic) {
    double delta = (update->ntimestep - time_origin) * dt;
    if (dynamic == VELOCITY) {
      x -= vx*delta;
      y -= vy*delta;
      z -= vz*delta;
    } else if (dynamic == WIGGLE) {
      double arg = omega_rotate * delta;
      double sine = sin(arg);
      x -= ax*sine;
      y -= ay*sine;
      z -= az*sine;
    } else if (dynamic == ROTATE) {
      double angle = -omega_rotate*delta;
      rotate(x,y,z,angle);
    }
  }

  return !(inside(x,y,z) ^ interior);
}

/* ----------------------------------------------------------------------
   generate list of contact points for interior or exterior regions
   if region is dynamic:
     change x by inverse of region change
     change contact point by region change
------------------------------------------------------------------------- */

int Region::surface(double x, double y, double z, double cutoff)
{
  int ncontact;
  double xnear[3],xhold[3];

  if (dynamic) {
    double delta = (update->ntimestep - time_origin) * dt;
    if (dynamic == VELOCITY) {
      x -= vx*delta;
      y -= vy*delta;
      z -= vz*delta;
    } else if (dynamic == WIGGLE) {
      double arg = omega_rotate * delta;
      double sine = sin(arg);
      x -= ax*sine;
      y -= ay*sine;
      z -= az*sine;
    } else if (dynamic == ROTATE) {
      xhold[0] = x;
      xhold[1] = y;
      xhold[2] = z;
      double angle = -omega_rotate*delta;
      rotate(x,y,z,angle);
    }
  }

  xnear[0] = x;
  xnear[1] = y;
  xnear[2] = z;

  if (interior) ncontact = surface_interior(xnear,cutoff);
  else ncontact = surface_exterior(xnear,cutoff);

  if (dynamic && ncontact) {
    double delta = (update->ntimestep - time_origin) * dt;
    if (dynamic == ROTATE) {
      for (int i = 0; i < ncontact; i++) {
	x -= contact[i].delx;
	y -= contact[i].dely;
	z -= contact[i].delz;
	double angle = omega_rotate*delta;
	rotate(x,y,z,angle);
	contact[i].delx = xhold[0] - x;
	contact[i].dely = xhold[1] - y;
	contact[i].delz = xhold[2] - z;
      }
    }
  }

  return ncontact;
}

/* ----------------------------------------------------------------------
   add a single contact at Nth location in contact array
   x = particle position
   xp,yp,zp = region surface point
------------------------------------------------------------------------- */

void Region::add_contact(int n, double *x, double xp, double yp, double zp)
{
  double delx = x[0] - xp;
  double dely = x[1] - yp;
  double delz = x[2] - zp;
  contact[n].r = sqrt(delx*delx + dely*dely + delz*delz);
  contact[n].delx = delx;
  contact[n].dely = dely;
  contact[n].delz = delz;
}

/* ----------------------------------------------------------------------
   rotate x,y,z by angle via right-hand rule around point and runit normal
   sign of angle determines whether rotating forward/backward in time
   return updated x,y,z
   P = point = vector = point of rotation
   R = vector = axis of rotation
   w = omega of rotation (from period)
   X0 = x,y,z = initial coord of atom
   R0 = runit = unit vector for R
   C = (X0 dot R0) R0 = projection of atom coord onto R
   D = X0 - P = vector from P to X0
   A = D - C = vector from R line to X0
   B = R0 cross A = vector perp to A in plane of rotation
   A,B define plane of circular rotation around R line
   x,y,z = P + C + A cos(w*dt) + B sin(w*dt)
------------------------------------------------------------------------- */

void Region::rotate(double &x, double &y, double &z, double angle)
{
  double a[3],b[3],c[3],d[3],disp[3];

  double sine = sin(angle);
  double cosine = cos(angle);
  double x0dotr = x*runit[0] + y*runit[1] + z*runit[2];
  c[0] = x0dotr * runit[0];
  c[1] = x0dotr * runit[1];
  c[2] = x0dotr * runit[2];
  d[0] = x - point[0];
  d[1] = y - point[1];
  d[2] = z - point[2];
  a[0] = d[0] - c[0];
  a[1] = d[1] - c[1];
  a[2] = d[2] - c[2];
  b[0] = runit[1]*a[2] - runit[2]*a[1];
  b[1] = runit[2]*a[0] - runit[0]*a[2];
  b[2] = runit[0]*a[1] - runit[1]*a[0];
  disp[0] = a[0]*cosine  + b[0]*sine;
  disp[1] = a[1]*cosine  + b[1]*sine;
  disp[2] = a[2]*cosine  + b[2]*sine;
  x = point[0] + c[0] + disp[0];
  y = point[1] + c[1] + disp[1];
  z = point[2] + c[2] + disp[2];
}

/* ---------------------------------------------------------------------- */

void Region::reset_random(int new_seed)
{
    if(comm->me == 0) fprintf(screen,"INFO: Resetting random generator for region %s\n",id);
    random->reset(new_seed);
}

/* ---------------------------------------------------------------------- */

void Region::generate_random(double *pos)
{
    if(!bboxflag) error->all("Impossible to generate random points on region with incomputable bounding box");
    do
    {
        pos[0] = extent_xlo + random->uniform() * (extent_xhi - extent_xlo);
        pos[1] = extent_ylo + random->uniform() * (extent_yhi - extent_ylo);
        pos[2] = extent_zlo + random->uniform() * (extent_zhi - extent_zlo);
    }
    while(!match(pos[0],pos[1],pos[2]));
}

/* ---------------------------------------------------------------------- */

// generates a random point that has a min distance from surface
void Region::generate_random_cut_away(double *pos,double cut)
{
    if(!bboxflag) error->all("Impossible to generate random points on region with incomputable bounding box");
    do
    {
        pos[0] = extent_xlo + random->uniform() * (extent_xhi - extent_xlo);
        pos[1] = extent_ylo + random->uniform() * (extent_yhi - extent_ylo);
        pos[2] = extent_zlo + random->uniform() * (extent_zhi - extent_zlo);
    }
    while(!match(pos[0],pos[1],pos[2]) || match_cut(pos,cut));
}

/* ---------------------------------------------------------------------- */

// generates a random point that has a min distance from surface
void Region::generate_random_within_cut(double *pos,double cut)
{
    if(!bboxflag) error->all("Impossible to generate random points on region with incomputable bounding box");
    do
    {
        pos[0] = extent_xlo + random->uniform() * (extent_xhi - extent_xlo);
        pos[1] = extent_ylo + random->uniform() * (extent_yhi - extent_ylo);
        pos[2] = extent_zlo + random->uniform() * (extent_zhi - extent_zlo);
    }
    while(!match_cut(pos,cut));
}

/* ---------------------------------------------------------------------- */

int Region::match_cut(double *pos,double cut)
{
  double a[3],b[3],c[3],d[3],x[3];
  vectorCopy3D(pos,x);

  if (dynamic) {
    double delta = (update->ntimestep - time_origin) * dt;
    if (dynamic == VELOCITY) {
      x[0] -= vx*delta;
      x[1] -= vy*delta;
      x[2] -= vz*delta;
    } else if (dynamic == WIGGLE) {
      double arg = omega_rotate * delta;
      double sine = sin(arg);
      x[0] -= ax*sine;
      x[1] -= ay*sine;
      x[2] -= az*sine;
    } else if (dynamic == ROTATE) {
      double angle = -omega_rotate*delta;
      rotate(x[0],x[1],x[2],angle);
    }
  }

  if(interior) return surface_interior(x,cut);
  else return surface_exterior(x,cut);
}

/* ---------------------------------------------------------------------- */

double Region::volume_mc(int n_test)
{
    // impossible to calculate volume if bbox non-existent
    if(!bboxflag) return 0;

    double pos[3],bbox_vol,region_vol;
    int n_in = 0;

    for(int i = 0; i < n_test; i++)
    {
        pos[0] = extent_xlo + random->uniform() * (extent_xhi - extent_xlo);
        pos[1] = extent_ylo + random->uniform() * (extent_yhi - extent_ylo);
        pos[2] = extent_zlo + random->uniform() * (extent_zhi - extent_zlo);
        if(!domain->is_in_subdomain(pos)) continue;
        if(match(pos[0],pos[1],pos[2])) n_in++;
    }

    MyMPI::My_MPI_Sum_Scalar(n_in,world);
    if(n_in == 0) error->all("Unable to calculate region volume - are you operating on a 2d region?");

    bbox_vol = (extent_xhi - extent_xlo) * (extent_yhi - extent_ylo) * (extent_zhi - extent_zlo);
    region_vol = static_cast<double>(n_in)/static_cast<double>(n_test) * bbox_vol;
    return region_vol;
}
