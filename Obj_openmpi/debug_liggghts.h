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

#ifndef LMP_DEBUG_LIGGGHTS_H
#define LMP_DEBUG_LIGGGHTS_H

#include "lammps.h"
#include "string.h"
#include "stdlib.h"
#include "style_fix.h"
#include "myvector.h"

namespace LAMMPS_NS {

inline void __debug__(LAMMPS* lmp)
{
    //fprintf(lmp->screen,"test");

    for(int i = 0; i < lmp->modify->nfix; i++)
    {
        if(strcmp(lmp->modify->fix[i]->style,"rigid/multisphere") == 0)
        {
            FixRigidMultisphere *fr = static_cast<FixRigidMultisphere*>(lmp->modify->fix[i]);
            int nb,np;
            double **omega = fr->get_dump_ref(nb,np,"omega");

            printVec3D(lmp->screen,"omega1",omega[1]);

        }

    }
}

}

#endif
