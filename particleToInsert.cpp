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

#include "particleToInsert.h"
#include "math.h"
#include "error.h"
#include "update.h"
#include "domain.h"
#include "atom.h"
#include "atom_vec.h"
#include "fix.h"
#include "myvector.h"
#include "modify.h"

using namespace LAMMPS_NS;

ParticleToInsert::ParticleToInsert(LAMMPS* lmp,int ns) : Pointers(lmp)
{
        groupbit = 0;

        nspheres = ns;

        x_ins = memory->create_2d_double_array(nspheres,3,"x_ins");
        radius_ins = new double[nspheres];
}

/* ---------------------------------------------------------------------- */

ParticleToInsert::~ParticleToInsert()
{
        memory->destroy_2d_double_array(x_ins);
        delete []radius_ins;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::insert()
{
    // perform the actual insertion
    // add particles, set coordinate and radius
    // set group mask to "all" plus fix groups

    int inserted = 0;

    if (domain->is_in_extended_subdomain(x_ins[0]))
    {
            inserted++;
            atom->avec->create_atom(atom_type,x_ins[0]);
            int m = atom->nlocal - 1;
            atom->mask[m] = 1 | groupbit;
            vectorCopy3D(v_ins,atom->v[m]);
            vectorCopy3D(omega_ins,atom->omega[m]);
            atom->radius[m] = radius_ins[0];
            atom->density[m] = density_ins;
            atom->rmass[m] = mass_ins;
            int nfix = modify->nfix;
            Fix **fix = modify->fix;
            for (int j = 0; j < nfix; j++)
               if (fix[j]->create_attribute) fix[j]->set_arrays(m);
    }
    return inserted;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::check_near_set_x_v_omega(double *x,double *v, double *omega, double *quat, double **xnear, int &nnear)
{
    // check sphere against all others in xnear
    // if no overlap add to xnear
    double del[3], rsq, radsum;

    vectorCopy3D(x,x_ins[0]);

    for(int i = 0; i < nnear; i++)
    {
        vectorSubtract3D(x_ins[0],xnear[i],del);
        rsq = vectorMag3DSquared(del);
        radsum = radius_ins[0] + xnear[i][3];

        // no success in overlap
        if (rsq <= radsum*radsum) return 0;
    }

    // no overlap with any other - success

    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    // add to xnear
    vectorCopy3D(x_ins[0],xnear[nnear]);
    xnear[nnear][3] = radius_ins[0];
    nnear++;

	return 1;
}

/* ---------------------------------------------------------------------- */

int ParticleToInsert::set_x_v_omega(double *x, double *v, double *omega, double *quat)
{

    vectorCopy3D(x,x_ins[0]);
    vectorCopy3D(v,v_ins);
    vectorCopy3D(omega,omega_ins);

    return 1;
}
