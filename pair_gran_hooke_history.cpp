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
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "pair_gran_hooke_history.h"
#include "atom.h"
#include "atom_vec.h"
#include "domain.h"
#include "force.h"
#include "update.h"
#include "modify.h"
#include "fix.h"
#include "fix_contact_history.h"
#include "comm.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "memory.h"
#include "error.h"
#include "fix_rigid.h"
#include "fix_property_global.h"
#include "mech_param_gran.h"
#include "compute_pair_gran_local.h"
#include "myvector.h"

using namespace LAMMPS_NS;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::PairGranHookeHistory(LAMMPS *lmp) : PairGran(lmp)
{
    //flag that we intend to use contact history
    history = 1;
    dnum = 3;

    Yeff = NULL;
    Geff = NULL;
    betaeff = NULL;
    veff = NULL;
    cohEnergyDens = NULL;
    coeffRestLog = NULL;
    coeffFrict = NULL;
    coeffRollFrict = NULL;

    charVelflag = 1;
}

/* ---------------------------------------------------------------------- */

PairGranHookeHistory::~PairGranHookeHistory()
{

}

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::history_args(char** args)
{
    //provide names and newtonflags for each history value
    //newtonflag = 0 means that the value
    args[0] = (char *) "shearx";
    args[1] = (char *) "1";
    args[2] = (char *) "sheary";
    args[3] = (char *) "1";
    args[4] = (char *) "shearz";
    args[5] = (char *) "1";
}

/* ---------------------------------------------------------------------- */
#define LMP_GRAN_DEFS_DEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_DEFINE

inline void PairGranHookeHistory::addCohesionForce(int &ip, int &jp,double &r, double &Fn_coh) 
{
    //r is the distance between the sphere's centeres
    double Acont = - M_PI/4 * ( (r-ri-rj)*(r+ri-rj)*(r-ri+rj)*(r+ri+rj) )/(r*r); //contact area of the two spheres
    
    Fn_coh=cohEnergyDens[itype][jtype]*Acont;
}

/* ---------------------------------------------------------------------- */

inline void PairGranHookeHistory::deriveContactModelParams(int &ip, int &jp,double &meff,double &deltan, double &kn, double &kt, double &gamman, double &gammat, double &xmu, double &rmu) 
{
    double reff=ri*rj/(ri+rj);
    kn = 16./15.*sqrt(reff)*(Yeff[itype][jtype])*pow(15.*meff*charVel*charVel/(16.*sqrt(reff)*Yeff[itype][jtype]),0.2);
    kt = kn;
    gamman=sqrt(4.*meff*kn/(1.+(M_PI/coeffRestLog[itype][jtype])*(M_PI/coeffRestLog[itype][jtype])));
    gammat=gamman;
    xmu=coeffFrict[itype][jtype];
    if(rollingflag)rmu=coeffRollFrict[itype][jtype];
    if (dampflag == 0) gammat = 0.0;

    // convert Kn and Kt from pressure units to force/distance^2
    kn /= force->nktv2p;
    kt /= force->nktv2p;

    return;
}

#define LMP_GRAN_DEFS_UNDEFINE
#include "pair_gran_defs.h"
#undef LMP_GRAN_DEFS_UNDEFINE

/* ---------------------------------------------------------------------- */

void PairGranHookeHistory::compute(int eflag, int vflag,int addflag)
{
  //calculated from the material properties 
  double kn,kt,gamman,gammat,xmu,rmu; 
  double Fn_coh;

  int i,j,ii,jj,inum,jnum,itype,jtype;
  double xtmp,ytmp,ztmp,delx,dely,delz,fx,fy,fz;
  double radi,radj,radsum,rsq,r,rinv,rsqinv,reff;
  double vr1,vr2,vr3,vnnr,vn1,vn2,vn3,vt1,vt2,vt3,wrmag;
  double wr1,wr2,wr3;
  double vtr1,vtr2,vtr3,vrel;
  double meff,damp,ccel,tor1,tor2,tor3,r_torque[3],r_torque_n[3];
  double fn,fs,fs1,fs2,fs3;
  double shrmag,rsht, cri, crj;
  int *ilist,*jlist,*numneigh,**firstneigh;
  int *touch,**firsttouch;
  double *shear,*allshear,**firstshear;

  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = vflag_fdotr = 0;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  double *mass = atom->mass;
  int *type = atom->type;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = listgranhistory->firstneigh;
  firstshear = listgranhistory->firstdouble;

  if (update->ntimestep > laststep) shearupdate = 1;
  else shearupdate = 0;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];
    touch = firsttouch[i];
    allshear = firstshear[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx*delx + dely*dely + delz*delz;
      radj = radius[j];
      radsum = radi + radj;

      if (rsq >= radsum*radsum) {

	// unset non-touching neighbors

        touch[jj] = 0;
        shear = &allshear[dnum*jj];
        shear[0] = 0.0;
        shear[1] = 0.0;
        shear[2] = 0.0;

      } else {
        r = sqrt(rsq);
        rinv = 1.0/r;
        rsqinv = 1.0/rsq;

        // relative translational velocity

        vr1 = v[i][0] - v[j][0];
        vr2 = v[i][1] - v[j][1];
        vr3 = v[i][2] - v[j][2];

        // normal component

        vnnr = vr1*delx + vr2*dely + vr3*delz;
        vn1 = delx*vnnr * rsqinv;
        vn2 = dely*vnnr * rsqinv;
        vn3 = delz*vnnr * rsqinv;

        // tangential component

        vt1 = vr1 - vn1;
        vt2 = vr2 - vn2;
        vt3 = vr3 - vn3;

        // relative rotational velocity

        double deltan=radsum-r;
        cri = radi-0.5*deltan;
        crj = radj-0.5*deltan;
        wr1 = (cri*omega[i][0] + crj*omega[j][0]) * rinv;
        wr2 = (cri*omega[i][1] + crj*omega[j][1]) * rinv;
        wr3 = (cri*omega[i][2] + crj*omega[j][2]) * rinv;

        // normal forces = Hookian contact + normal velocity damping

        double mi,mj;
        if (rmass) {
          mi = rmass[i];
          mj = rmass[j];
        } else {
          itype = type[i];
          jtype = type[j];
          mi = mass[itype];
          mj = mass[jtype];
        }
        if (fr)
        {
           if(fr->body[i] >= 0) mi = fr->masstotal[fr->body[i]];  
           if(fr->body[j] >= 0) mj = fr->masstotal[fr->body[j]];  
        }

        meff = mi*mj/(mi+mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        deriveContactModelParams(i,j,meff,deltan,kn,kt,gamman,gammat,xmu,rmu);	 //modified C.K

        damp = gamman*vnnr*rsqinv;  
        ccel = kn*(radsum-r)*rinv - damp;
        
        if (cohesionflag) { 
            addCohesionForce(i,j,r,Fn_coh);
            ccel-=Fn_coh*rinv;
        }

        // relative velocities

        vtr1 = vt1 - (delz*wr2-dely*wr3);
        vtr2 = vt2 - (delx*wr3-delz*wr1);
        vtr3 = vt3 - (dely*wr1-delx*wr2);
        vrel = vtr1*vtr1 + vtr2*vtr2 + vtr3*vtr3;
        vrel = sqrt(vrel);

        // shear history effects

        touch[jj] = 1;

        shear = &allshear[dnum*jj];

        if (shearupdate && addflag)
        {
            shear[0] += vtr1*dt;
            shear[1] += vtr2*dt;
            shear[2] += vtr3*dt;

            // rotate shear displacements

            rsht = shear[0]*delx + shear[1]*dely + shear[2]*delz;
            rsht *= rsqinv;
            shear[0] -= rsht*delx;
            shear[1] -= rsht*dely;
            shear[2] -= rsht*delz;
        }

        shrmag = sqrt(shear[0]*shear[0] + shear[1]*shear[1] +  shear[2]*shear[2]);

        // tangential forces = shear + tangential velocity damping

        fs1 = - (kt*shear[0]);
        fs2 = - (kt*shear[1]);
        fs3 = - (kt*shear[2]);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1*fs1 + fs2*fs2 + fs3*fs3);
        fn = xmu * fabs(ccel*r);

        // energy loss from sliding or damping
        if (fs > fn) {
            if (shrmag != 0.0) {
                fs1 *= fn/fs;
                fs2 *= fn/fs;
                fs3 *= fn/fs;
                shear[0] = -fs1/kt;
                shear[1] = -fs2/kt;
                shear[2] = -fs3/kt;
            }
            else fs1 = fs2 = fs3 = 0.0;
        }
        else
        {
            fs1 -= (gammat*vtr1);
            fs2 -= (gammat*vtr2);
            fs3 -= (gammat*vtr3);
        }

        // forces & torques

        fx = delx*ccel + fs1;
        fy = dely*ccel + fs2;
        fz = delz*ccel + fs3;

        tor1 = rinv * (dely*fs3 - delz*fs2);
        tor2 = rinv * (delz*fs1 - delx*fs3);
        tor3 = rinv * (delx*fs2 - dely*fs1);

        // add rolling friction torque
        vectorZeroize3D(r_torque);
        if(rollingflag)
        {
            reff=radi*radj/(radi+radj);
            wrmag = sqrt(wr1*wr1+wr2*wr2+wr3*wr3);
            if (wrmag > 0.)
            {
                r_torque[0] = rmu*kn*deltan*wr1/wrmag*reff;
                r_torque[1] = rmu*kn*deltan*wr2/wrmag*reff;
                r_torque[2] = rmu*kn*deltan*wr3/wrmag*reff;
            }
            // remove normal (torsion) part of torque
            r_torque_n[0] = r_torque[0] * delx * rinv;
            r_torque_n[1] = r_torque[1] * dely * rinv;
            r_torque_n[2] = r_torque[2] * delz * rinv;
            vectorSubtract3D(r_torque,r_torque_n,r_torque);
        }

        if(addflag)
        {
            f[i][0] += fx;
            f[i][1] += fy;
            f[i][2] += fz;
            torque[i][0] -= cri*tor1 + r_torque[0];
            torque[i][1] -= cri*tor2 + r_torque[1];
            torque[i][2] -= cri*tor3 + r_torque[2];
        }

        if (j < nlocal && addflag) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;
          torque[j][0] -= crj*tor1 + r_torque[0];
          torque[j][1] -= crj*tor2 + r_torque[1];
          torque[j][2] -= crj*tor3 + r_torque[2];
        }

        if(cpl && !addflag) cpl->add_pair(i,j,fx,fy,fz,tor1,tor2,tor3,shear);

        if (evflag) ev_tally_xyz(i,j,nlocal,0,0.0,0.0,fx,fy,fz,delx,dely,delz);
      }
    }
  }

  laststep = update->ntimestep;
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHookeHistory::settings(int narg, char **arg) 
{
  if (narg != 2) error->all("Illegal pair_style command");

  dampflag = force->inumeric(arg[0]) & 1;
  rollingflag = force->inumeric(arg[0]) & 2;
  cohesionflag = force->inumeric(arg[1]);

  if (dampflag < 0 || dampflag > 3 || cohesionflag < 0 || cohesionflag > 1)
    error->all("Illegal pair_style command");

  if(cohesionflag && domain->dimension!=3) error->all("Cohesion model valid for 3d simulations only");
}

/* ----------------------------------------------------------------------
   init specific to this granular substyle
------------------------------------------------------------------------- */

void PairGranHookeHistory::init_substyle()
{
  int max_type = mpg->max_type();

  allocate_properties(max_type);

  //Get pointer to the fixes that have the material properties
  
  Y1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("youngsModulus","property/global","peratomtype",max_type,0));
  v1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("poissonsRatio","property/global","peratomtype",max_type,0));

  coeffRest1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRestitution","property/global","peratomtypepair",max_type,max_type));
  coeffFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientFriction","property/global","peratomtypepair",max_type,max_type));

  if(rollingflag)
    coeffRollFrict1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("coefficientRollingFriction","property/global","peratomtypepair",max_type,max_type));
  if(cohesionflag)
    cohEnergyDens1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("cohesionEnergyDensity","property/global","peratomtypepair",max_type,max_type));

  if(charVelflag) charVel1=static_cast<FixPropertyGlobal*>(modify->find_fix_property("characteristicVelocity","property/global","scalar",0,0));

  //pre-calculate parameters for possible contact material combinations
  for(int i=1;i< max_type+1; i++)
  {
      for(int j=1;j<max_type+1;j++)
      {
          double Yi=Y1->compute_vector(i-1);
          double Yj=Y1->compute_vector(j-1);
          double vi=v1->compute_vector(i-1);
          double vj=v1->compute_vector(j-1);

          Yeff[i][j] = 1./((1.-pow(vi,2.))/Yi+(1.-pow(vj,2.))/Yj);
          Geff[i][j] = 1./(2.*(2.-vi)*(1.+vi)/Yi+2.*(2.-vj)*(1.+vj)/Yj);

          coeffRestLog[i][j] = log(coeffRest1->compute_array(i-1,j-1));

          betaeff[i][j] =coeffRestLog[i][j] /sqrt(pow(coeffRestLog[i][j],2.)+pow(M_PI,2.));

          coeffFrict[i][j] = coeffFrict1->compute_array(i-1,j-1);
          if(rollingflag) coeffRollFrict[i][j] = coeffRollFrict1->compute_array(i-1,j-1);

          if(cohesionflag) cohEnergyDens[i][j] = cohEnergyDens1->compute_array(i-1,j-1);
          //omitting veff here

      }
  }

  if(charVelflag) charVel = charVel1->compute_scalar();
}

/* ----------------------------------------------------------------------
  allocate per-type and per-type pair properties
------------------------------------------------------------------------- */

void PairGranHookeHistory::allocate_properties(int size)
{
    memory->destroy_2d_double_array(Yeff);
    memory->destroy_2d_double_array(Geff);
    memory->destroy_2d_double_array(betaeff);
    memory->destroy_2d_double_array(veff);
    memory->destroy_2d_double_array(cohEnergyDens);
    memory->destroy_2d_double_array(coeffRestLog);
    memory->destroy_2d_double_array(coeffFrict);
    memory->destroy_2d_double_array(coeffRollFrict);
    Yeff = memory->create_2d_double_array(size+1,size+1,"Yeff");
    Geff = memory->create_2d_double_array(size+1,size+1,"Geff");
    betaeff = memory->create_2d_double_array(size+1,size+1,"betaeff");
    veff = memory->create_2d_double_array(size+1,size+1,"veff");
    cohEnergyDens = memory->create_2d_double_array(size+1,size+1,"cohEnergyDens");
    coeffRestLog = memory->create_2d_double_array(size+1,size+1,"coeffRestLog");
    coeffFrict = memory->create_2d_double_array(size+1,size+1,"coeffFrict");
    coeffRollFrict = memory->create_2d_double_array(size+1,size+1,"coeffRollFrict");
}

/* ----------------------------------------------------------------------
  proc 0 writes to restart file 
------------------------------------------------------------------------- */

void PairGranHookeHistory::write_restart_settings(FILE *fp) 
{
  int writeflag = dampflag + rollingflag * 2;
  fwrite(&writeflag,sizeof(int),1,fp);
  fwrite(&cohesionflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
  proc 0 reads from restart file, bcasts 
------------------------------------------------------------------------- */

void PairGranHookeHistory::read_restart_settings(FILE *fp) 
{
  if (comm->me == 0) {
    int readflag;
    fread(&readflag,sizeof(int),1,fp);
    fread(&cohesionflag,sizeof(int),1,fp);
    dampflag = readflag & 1;
    rollingflag = readflag & 2;
  }
  MPI_Bcast(&dampflag,1,MPI_INT,0,world);
  MPI_Bcast(&cohesionflag,1,MPI_INT,0,world);
  MPI_Bcast(&rollingflag,1,MPI_INT,0,world);
}
