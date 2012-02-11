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
#include "fix_insert_stream.h"
#include "fix_mesh_gran.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "update.h"
#include "comm.h"
#include "modify.h"
#include "myvector.h"
#include "domain.h"
#include "random_park.h"
#include "memory.h"
#include "error.h"
#include "fix_property_atom.h"
#include "fix_particledistribution_discrete.h"
#include "fix_template_sphere.h"
#include "particleToInsert.h"

enum{FACE_NONE,FACE_MESH,FACE_CIRCLE};

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixInsertStream::FixInsertStream(LAMMPS *lmp, int narg, char **arg) :
  FixInsert(lmp, narg, arg)
{
  // set defaults first, then parse args
  init_defaults();

  bool hasargs = true;
  while(iarg < narg && hasargs)
  {
    hasargs = false;
    if (strcmp(arg[iarg],"insertion_face") == 0) {
      
      if (iarg+2 > narg) error->all("Illegal fix insert/stream command");
      int f_i = modify->find_fix(arg[iarg+1]);
      if (f_i == -1) error->all("Could not find fix mesh/gran id you provided for the fix insert/stream command");
      if (strncmp(modify->fix[f_i]->style,"mesh/gran",8)) error->all("Fix insert/stream: The fix belonging to the id you provided is not of type mesh/gran");
      ins_face = static_cast<FixMeshGran*>(modify->fix[f_i]);
      face_style = FACE_MESH;
      iarg += 2;
      hasargs = true;
    }else if (strcmp(arg[iarg],"extrude_length") == 0) {
      if (iarg+2 > narg) error->all("Illegal fix insert/stream command");
      extrude_length = atof(arg[iarg+1]);
      if(extrude_length < 0. ) error->all("Illegal fix insert/stream command, invalid extrude_length");
      iarg += 2;
      hasargs = true;
    }else error->all("Illegal fix insert/stream command, unknown keyword");
  }

  fix_release = NULL;

  fix_rm = NULL;

  releasevar_name = new char[10+strlen(id)];

  nevery = 1;
}

/* ---------------------------------------------------------------------- */

FixInsertStream::~FixInsertStream()
{
    delete []releasevar_name;

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::post_create()
{
  FixInsert::post_create();

  if(!fix_release)
  {
        char* fixarg[13];

        fixarg[0] = new char[10+strlen(id)];
        sprintf(fixarg[0],"release_%s",id);
        
        fixarg[1]="all";
        fixarg[2]="property/atom";
        fixarg[3] = new char[10+strlen(id)];
        sprintf(releasevar_name,"release_%s",id);
        strcpy(fixarg[3],releasevar_name);
        fixarg[4]="vector"; 
        fixarg[5]="yes";    
        fixarg[6]="yes";    
        fixarg[7]="no";    
        fixarg[8]="0.";
        fixarg[9]="0.";
        fixarg[10]="0.";
        fixarg[11]="0.";
        fixarg[12]="0.";
        fix_release = modify->add_fix_property_atom(13,fixarg);

        delete []fixarg[0];
        delete []fixarg[3];
  }
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init_defaults()
{
    face_style = FACE_NONE;
    extrude_length = 0.;
    insert_every = -1;
}

/* ----------------------------------------------------------------------
   calculate ninsert, insert_every, ninsert_per, massinsert, flowrates etc
   also perform error checks
------------------------------------------------------------------------- */

void FixInsertStream::calc_insertion_properties()
{
    double dt,dot,extrude_vec[3],t1[3],t2[3];

    // error check on insertion face
    if(face_style == FACE_NONE) error->all("Illegal fix insert/stream command, must define an insertion face");
    if(all_in_flag == 1) error->all("Illegal fix insert/stream command, all_in not available, must be set to false");

    // check properties of insertion face
    if(face_style == FACE_MESH)
    {
        // check if face planar
        if(!ins_face->STLdata->is_planar()) error->all("Fix insert/stream command requires a planar face for insertion");

        // get normal vector of face 0
        ins_face->STLdata->normal_vec(normalvec,0);

        // flip normal vector so dot product with v_insert is >0
        dot = vectorDot3D(v_insert,normalvec);
        if(dot < 0) vectorScalarMult3D(normalvec,-1.);

        // calc v normal
        dot = vectorDot3D(v_insert,normalvec);
        vectorCopy3D(normalvec,v_normal);
        vectorScalarMult3D(v_normal,dot);

        // error check on v normal
        if(vectorMag3D(v_normal) < 1.e-3) error->all("Illegal fix insert/stream command, insertion velocity projected on face normal is < 1e-3");

        // get reference point = com
        vectorCopy3D(ins_face->p_ref,p_ref);
    }
    else error->all("FixInsertStream::calc_insertion_properties(): Implementation missing");

    // error check on insertion velocity
    if(vectorMag3D(v_insert) < 1e-5) error->all("Illegal fix insert/stream command, insertion velocity too low");

    if(insert_every == -1 && extrude_length == 0.) error->all("Illegal fix insert/stream command, must define either 'insert_every' or 'extrude_length'");
    if(insert_every > -1 && extrude_length > 0.) error->all("Illegal fix insert/stream command, must not provide both 'insert_every' and 'extrude_length'");

    dt = update->dt;

    // if extrude_length given, calculate insert_every
    if(insert_every == -1)
    {
        if(extrude_length < 3.*max_r_bound()) error->all("Illegal fix insert/stream command, 'extrude_length' is too small");
        insert_every = static_cast<int>(extrude_length/(dt*vectorMag3D(v_normal)));
        
        if(insert_every == 0) error->all("Illegal fix insert/stream command, insertion velocity too high or extrude_length too low");
    }
    // if insert_every given, calculate extrude_length
    else
    {
        if(insert_every < 1) error->all("Illegal fix insert/stream command, 'insert_every' must be > 0");
        extrude_length = static_cast<double>(insert_every) * dt * vectorMag3D(v_normal);
        if(extrude_length < 3.*max_r_bound()) error->all("Illegal fix insert/stream command, 'insert_every' or 'vel' is too small");
    }

    // ninsert
    // if ninsert not defined directly, calculate it
    if(ninsert == 0)
    {
        if(massinsert > 0.) ninsert = static_cast<int>(massinsert / fix_distribution->mass_expect());
        else error->all("Illegal fix insert/stream command, must define either 'nparticles' or 'mass'");
    }

    // flow rate, ninsert_per
    if(nflowrate == 0.)
    {
        if(massflowrate == 0.) error->all("Illegal fix insert/stream command, must define either 'massrate' or 'particlerate'");
        nflowrate = massflowrate / fix_distribution->mass_expect();
    }
    else massflowrate = nflowrate * fix_distribution->mass_expect();

    ninsert_per = nflowrate*(static_cast<double>(insert_every)*dt);

    massinsert = static_cast<double>(ninsert) * fix_distribution->mass_expect();

    // calculate bounding box of extruded face
    if(face_style == FACE_MESH)
    {
        // get bounding box for face
        ins_face->STLdata->getMeshBoundBox(ins_vol_xmin,ins_vol_xmax);

        // get bounding box for extruded face - store in t1,t2
        vectorScalarMult3D(normalvec,-extrude_length,extrude_vec);
        vectorAdd3D(ins_vol_xmin,extrude_vec,t1);
        vectorAdd3D(ins_vol_xmax,extrude_vec,t2);

        // take min and max
        vectorComponentMin3D(ins_vol_xmin,t1,ins_vol_xmin);
        vectorComponentMax3D(ins_vol_xmax,t2,ins_vol_xmax);

    }
    else error->all("Missing implementation in calc_insertion_properties()");

}

/* ---------------------------------------------------------------------- */

int FixInsertStream::setmask()
{
    int mask = FixInsert::setmask();
    mask |= END_OF_STEP;
    return mask;
}

/* ---------------------------------------------------------------------- */

void FixInsertStream::init()
{
    
    FixInsert::init();

    for(int ifix = 0; ifix < modify->nfix; ifix++)
    {
        Fix *fix = modify->fix[ifix];
        if(strncmp(fix->style,"rigid/multisphere",17) == 0)
            fix_rm = (FixRigidMultisphere*) fix;
    }

    fix_release = static_cast<FixPropertyAtom*>(modify->find_fix_property(releasevar_name,"property/atom","vector",5,0));

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::pre_insert()
{
    if(!domain->is_in_domain(ins_vol_xmin) || !domain->is_in_domain(ins_vol_xmax))
        error->warning("Fix insert/stream: Extruded insertion face extends outside domain, may not insert all particles correctly");
}

/* ---------------------------------------------------------------------- */

inline int FixInsertStream::overlap(int i)
{
    double pos_rel[3], pos_projected[3], t[3];
    double **x = atom->x;

    vectorSubtract3D(x[i],p_ref,pos_rel);
    double dist_normal = vectorDot3D(pos_rel,normalvec);

    // on wrong side of extrusion
    if(dist_normal > 0) return 0;

    // on right side of extrusion, but too far away
    if(dist_normal < -extrude_length) return 0;

    // on right side of extrusion, within extrude_length
    // check if projection is on face or not

    vectorScalarMult3D(normalvec,dist_normal,t);
    vectorAdd3D(x[i],t,pos_projected);

    return ins_face->STLdata->is_on_surface(pos_projected);
}

/* ----------------------------------------------------------------------
   generate random positions on insertion face
   extrude by random length in negative face normal direction
     currently only implemented for all_in_flag = 0
     since would be tedious to check/implement otherwise
------------------------------------------------------------------------- */

inline void FixInsertStream::generate_random(double *pos, double rad)
{
    // generate random position on the mesh
    ins_face->STLdata->generate_random(pos,random);

    // extrude the position
    
    double r = -1.*(random->uniform()*(extrude_length-2.*rad) + rad);

    double ext[3];
    vectorScalarMult3D(normalvec,r,ext);
    vectorAdd3D(pos,ext,pos);
}

/* ----------------------------------------------------------------------
   generate random positions within extruded face
   perform overlap check via xnear if requested
   returns # bodies and # spheres that could actually be inserted
------------------------------------------------------------------------- */

void FixInsertStream::x_v_omega(int ninsert_this,int &ninserted_this, int &ninserted_spheres_this, double &mass_inserted_this)
{
    
    ninserted_this = ninserted_spheres_this = 0;
    mass_inserted_this = 0.;

    int nins;
    double pos[3];
    ParticleToInsert *pti;

    double omega_tmp[] = {0.,0.,0.};

    // no overlap check
    // insert with v_normal, no omega
    if(!check_ol_flag)
    {
        for(int itotal = 0; itotal < ninsert_this; itotal++)
        {
            pti = fix_distribution->pti_list[ninserted_this];
            double rad_to_insert = pti->r_bound_ins;
            generate_random(pos,rad_to_insert);

            // could ramdonize vel, omega, quat here

            nins = pti->set_x_v_omega(pos,v_normal,omega_tmp,quat_insert);

            ninserted_spheres_this += nins;
            mass_inserted_this += pti->mass_ins;
            ninserted_this++;
        }
    }
    // overlap check
    // account for maxattempt
    // pti checks against xnear and adds self contributions
    else
    {
        int ntry = 0;
        int maxtry = ninsert_this * maxattempt;

        while(ntry < maxtry && ninserted_this < ninsert_this)
        {
            pti = fix_distribution->pti_list[ninserted_this];
            double rad_to_insert = pti->r_bound_ins;

            nins = 0;
            while(nins == 0 && ntry < maxtry)
            {
                generate_random(pos,rad_to_insert);

                // could ramdonize vel, omega, quat here

                ntry++;
                nins = pti->check_near_set_x_v_omega(pos,v_normal,omega_tmp,quat_insert,xnear,nspheres_near);
            }

            if(nins > 0)
            {
                ninserted_spheres_this += nins;
                mass_inserted_this += pti->mass_ins;
                ninserted_this++;
            }
        }
    }

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::finalize_insertion(int nins)
{
    // nins particles have been inserted on this proc, set initial position, insertion step and release step according to pos

    int n_steps;
    int step = update->ntimestep;
    int ilo = atom->nlocal - nins;
    int ihi = atom->nlocal;

    double pos_rel[3], dist_normal;
    double **x = atom->x;
    double dt = update->dt;

    double **release_data = fix_release->array_atom;

    for(int i = ilo; i < ihi; i++)
    {
        
        vectorSubtract3D(p_ref,x[i],pos_rel);
        dist_normal = vectorDot3D(pos_rel,normalvec);
        n_steps = static_cast<int>(dist_normal/(vectorMag3D(v_normal)*dt));

        // first 3 values is original position to integrate
        vectorCopy3D(x[i],release_data[i]);

        // 4th value is insertion step
        release_data[i][3] = static_cast<double>(step);

        // 5th value is step to release
        release_data[i][4] = static_cast<double>(step + n_steps);

    }

}

/* ---------------------------------------------------------------------- */

void FixInsertStream::end_of_step()
{
    int r_step, i_step;

    int step = update->ntimestep;
    int nlocal = atom->nlocal;
    double **release_data = fix_release->array_atom;
    double time_elapsed, dist_elapsed[3];
    double dt = update->dt;

    double **x = atom->x;
    double **v = atom->v;
    double **f = atom->f;
    double **omega = atom->omega;
    double **torque = atom->torque;
    int *mask = atom->mask;

    for(int i = 0; i < nlocal; i++)
    {
        if (mask[i] & groupbit)
        {
            i_step = static_cast<int>(release_data[i][3]);
            r_step = static_cast<int>(release_data[i][4]);

            if(step > r_step) continue;
            else if (r_step == step)
            {
                // set inital conditions
                vectorCopy3D(v_insert,v[i]);
                vectorCopy3D(omega_insert,omega[i]);

                if(fix_rm) error->all("must set params for frm");

            }
            // step < r_step, only true for inserted particles
            //   b/c r_step is 0 for all other particles
            // integrate with constant vel
            else
            {
                if(fix_rm) error->all("must do this for frm");

                time_elapsed = (step - i_step) * dt;

                // every particle moves with v_normal
                vectorScalarMult3D(v_normal,time_elapsed,dist_elapsed);
                double *x_ins = release_data[i];

                // set x,v,omega
                vectorAdd3D(x_ins,dist_elapsed,x[i]);
                vectorCopy3D(v_normal,v[i]);
                vectorZeroize3D(omega[i]);

                // zero out force, torque
                vectorZeroize3D(f[i]);
                vectorZeroize3D(torque[i]);
            }
        }
    }

}
