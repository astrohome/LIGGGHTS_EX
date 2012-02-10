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

// contributing author for rotation option: Evan Smuts (U Cape Town)

#include "stl_tri.h"
#include "lammps.h"
#include "memory.h"
#include "error.h"
#include "math.h"
#include "domain.h"
#include "myvector.h"
#include "math_extra_liggghts.h"
#include "random_park.h"

#define BIG 1.0e20
#define SMALL 1.e-10
#define DELTATRI 5000 

using namespace LAMMPS_NS;

STLtri::STLtri(LAMMPS *lmp): Pointers(lmp)
{
    node=NULL;
    node_lastRe=NULL;
    v_node=NULL;
    f_tri=NULL;
    fn_fshear=NULL;
    facenormal=NULL;
    cK=NULL;
    ogK=NULL;
    ogKlen=NULL;
    oKO=NULL;
    rK=NULL;
    x=NULL;
    v=NULL;
    xoriginal=NULL;
    f=NULL;
    rmass=NULL;
    Area=NULL;
    wear=NULL;
    wear_step=NULL;
    neighfaces=NULL;
    contactInactive=NULL;

    nTri=0;
    nTriMax=0;
    xvf_len=0;
    xvf_lenMax=0;
    grow_arrays();
    before_rebuild(); 
    movingMesh=0;
    conveyor=0;
    rotation=0;
    skinSafetyFactor=0.;

    //bitmask settings
    EDGE_INACTIVE = new int[3];
    CORNER_INACTIVE = new int[3];
      EDGE_INACTIVE[0] = 1;    EDGE_INACTIVE[1] = 2;    EDGE_INACTIVE[2] = 4;
    CORNER_INACTIVE[0] = 8;  CORNER_INACTIVE[1] = 16; CORNER_INACTIVE[2] = 32;
}

/* ---------------------------------------------------------------------- */

STLtri::~STLtri()
{
   memory->destroy_3d_double_array(node);
   memory->destroy_2d_double_array(facenormal);
   memory->destroy_2d_double_array(f_tri);
   memory->destroy_2d_double_array(fn_fshear);
   memory->destroy_2d_double_array(cK);
   memory->destroy_3d_double_array(ogK);
   memory->destroy_2d_double_array(ogKlen);
   memory->destroy_3d_double_array(oKO);
   memory->sfree(rK);
   memory->destroy_2d_double_array(v);
   memory->destroy_2d_double_array(xoriginal);
   memory->destroy_2d_double_array(f);
   memory->sfree(rmass);
   memory->sfree(Area);
   memory->sfree(wear);
   memory->sfree(wear_step);
   memory->destroy_2d_int_array(neighfaces);
   memory->sfree(contactInactive);
   
   if(conveyor || rotation) memory->destroy_3d_double_array(v_node);

   delete []EDGE_INACTIVE;delete []CORNER_INACTIVE;
}

/* ---------------------------------------------------------------------- */

void STLtri::grow_arrays()
{
    nTriMax+=DELTATRI;
    xvf_lenMax+=DELTATRI*VECPERTRI;

    node=(double***)(memory->grow_3d_double_array(node,nTriMax+1, 3 , 3, "stl_tri_node"));
    node_lastRe=(double***)(memory->grow_3d_double_array(node_lastRe,nTriMax+1, 3 , 3, "stl_tri_node_lastRe"));
    facenormal=(double**)(memory->grow_2d_double_array(facenormal,nTriMax+1, 3, "stl_tri_facenormal"));
    
    cK=(double**)(memory->grow_2d_double_array(cK, nTriMax+1, 3, "stl_tri_cK"));
    ogK=(double***)(memory->grow_3d_double_array(ogK, nTriMax+1, 3, 3, "stl_tri_ogK"));
    ogKlen=(double**)(memory->grow_2d_double_array(ogKlen, nTriMax+1, 3, "stl_tri_ogKlen"));
    oKO=(double***)(memory->grow_3d_double_array(oKO, nTriMax+1, 3, 3, "stl_tri_oKO"));
    rK=(double*)(memory->srealloc(rK, (nTriMax+1)*sizeof(double), "stl_tri_rK"));
    Area=(double*)(memory->srealloc(Area, (nTriMax+1)*sizeof(double), "stl_tri_Area"));
    wear=(double*)(memory->srealloc(wear, (nTriMax+1)*sizeof(double), "stl_tri_wear"));
    wear_step=(double*)(memory->srealloc(wear_step, (nTriMax+1)*sizeof(double), "stl_tri_wear_step"));
    f_tri=(double**)(memory->grow_2d_double_array(f_tri,(nTriMax+1),3, "stl_tri_f_tri"));
    fn_fshear=(double**)(memory->grow_2d_double_array(fn_fshear,(nTriMax+1),3,"stl_tri_fn_fshear"));
    neighfaces=(int**)(memory->grow_2d_int_array(neighfaces,(nTriMax+1), 3, "stl_tri_neighfaces"));
    contactInactive=(int*)(memory->srealloc(contactInactive, (nTriMax+1)*sizeof(int), "stl_tri_contactActive"));

    for (int i=0;i<(nTriMax+1);i++) {
        f_tri[i][0]=0.;f_tri[i][2]=0.;f_tri[i][2]=0.;
        fn_fshear[i][0]=0.;fn_fshear[i][1]=0.;fn_fshear[i][2]=0.;
        wear[i] = 0.;

        neighfaces[i][0]=-1;neighfaces[i][1]=-1;neighfaces[i][2]=-1;
        contactInactive[i]=0; //init with all contacts active
    }
}

/* ---------------------------------------------------------------------- */
//initialize conveyor model
void STLtri::initConveyor()
{
    conveyor=1;
    if(v_node==NULL) v_node=(double***)(memory->grow_3d_double_array(v_node,(nTriMax+1), 3 , 3, "stl_tri_v_node"));
    double vtri[3];
    double tmp[3];
    double scp;
    double conv_vel_mag=vectorMag3D(conv_vel);

    for (int i=0;i<nTri;i++)
    {

        scp=vectorDot3D(conv_vel,facenormal[i]);
        vectorScalarMult3D(facenormal[i],scp,tmp);
        for(int j=0;j<3;j++)
        {
            
            vectorSubtract3D(conv_vel,tmp,v_node[i][j]);
            if(vectorMag3D(v_node[i][j])>0.)
            {
                vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
                vectorScalarMult3D(v_node[i][j],conv_vel_mag);
            }
           
        }
    }
}

/* ---------------------------------------------------------------------- */

//initialize rotation model
void STLtri::initRotation(double epsilon)  	//added by evan
{
    rotation=1;
    if(v_node==NULL) v_node=(double***)(memory->grow_3d_double_array(v_node,(nTriMax+1), 3 , 3, "stl_tri_v_node"));
    double tmp[3];
    double scp;
    double unitAxis[3];
    double tangComp[3];
    double Utang[3];
    double surfaceV[3];

    double magAxis = vectorMag3D(rot_axis);
    vectorScalarDiv3D(rot_axis,magAxis,unitAxis); 		//calculate unit vector of rotation axis

    for (int i=0;i<nTri;i++)	//number of STL faces
    {
        for(int j=0;j<3;j++) 	//number of nodes per face (3)
        {
	    vectorSubtract3D(node[i][j],rot_origin,surfaceV); 	//position of node - origin of rotation (to get lever arm)
    	    vectorCross3D(surfaceV,unitAxis,tangComp);	      	//lever arm X rotational axis = tangential component
    	    vectorScalarMult3D(tangComp,-rot_omega,Utang);     	//multiplying by omega scales the tangential component to give tangential velocity
	    if(vectorMag3D(Utang)<epsilon) error->all("Rotation velocity too low"); //EPSILON is wall velocity, not rotational omega
            scp = vectorDot3D(Utang, facenormal[i]);
            vectorScalarMult3D(facenormal[i],scp,tmp);
            vectorSubtract3D(Utang,tmp,v_node[i][j]);	      //removes components normal to wall
    	    double magUtang = vectorMag3D(Utang);
            if(vectorMag3D(v_node[i][j])>0.)
            {
               vectorScalarDiv3D(v_node[i][j],vectorMag3D(v_node[i][j]));
               vectorScalarMult3D(v_node[i][j],magUtang);
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

//copy values of x to xoriginal
void STLtri::initMove(double skinSafety)
{
    if(conveyor || rotation) error->all("Can not use moving mesh feature together with conveyor feature");
    movingMesh=1;
    skinSafetyFactor=skinSafety;

    xvf_len=nTri*VECPERTRI;

    if(v_node==NULL) v_node=(double***)(memory->grow_3d_double_array(v_node,(nTriMax+1), 3 , 3, "stl_tri_v_node"));
    if(x==NULL) x = (double**)(memory->srealloc(x,xvf_lenMax*sizeof(double *),"stl_tri:x"));
    if(v==NULL) v = (double**)(memory->grow_2d_double_array(v, xvf_lenMax, 3, "stl_tri:v"));

    //zero node velocity
    for (int i=0;i<xvf_lenMax;i++){
        for (int j=0;j<3;j++) v[i][j]=0.;
    }

    int m=0;
    for (int i=0;i<nTriMax;i++)
    {
        
        for (int j=0;j<3;j++) {
            v_node[i][j]=v[m];
            x[m++]=node[i][j];
        }

        x[m++]=facenormal[i];
        x[m++]=cK[i];
        for (int j=0;j<3;j++) x[m++]=ogK[i][j];
        for (int j=0;j<3;j++) x[m++]=oKO[i][j];
    }
    if(xoriginal==NULL) xoriginal = (double**)(memory->grow_2d_double_array(xoriginal, xvf_lenMax, 3, "stl_tri:xoriginal"));
    if(f==NULL) f = (double**)(memory->grow_2d_double_array(f, xvf_lenMax, 3, "stl_tri:f"));
    if(rmass==NULL) rmass=(double*)(memory->srealloc(rmass, xvf_lenMax*sizeof(double), "stl_import_rmass"));

    vecToPoint();
    for (int i=0;i<xvf_len;i++)
    {
       for(int j=0;j<3;j++) xoriginal[i][j]=x[i][j];
       rmass[i]=1.;
    }
    pointToVec();

}

/* ---------------------------------------------------------------------- */

//save nodes positions
void STLtri::before_rebuild()
{
    for (int i=0;i<nTri;i++)
    {
        for(int j=0;j<3;j++)
        {
            vectorCopy3D(node[i][j],node_lastRe[i][j]);
        }
    }
}

/* ---------------------------------------------------------------------- */

void STLtri::vecToPoint()
{
    for (int i=0;i<nTri;i++)
    {
        
        for(int j=0;j<3;j++) {
            facenormal[i][j]+=cK[i][j];
        }

        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                oKO[i][j][k]+=node[i][j][k];
                ogK[i][j][k]+=cK[i][k];
            }
        }
    }
}

/* ---------------------------------------------------------------------- */

void STLtri::pointToVec()
{
    for (int i=0;i<nTri;i++)
    {
        
        for(int j=0;j<3;j++) facenormal[i][j]-=cK[i][j];

        vectorScalarDiv3D(facenormal[i],vectorMag3D(facenormal[i]));

        for(int j=0;j<3;j++){
            for(int k=0;k<3;k++){
                oKO[i][j][k]-=node[i][j][k];
                ogK[i][j][k]-=cK[i][k];
            }
            
            vectorScalarDiv3D(oKO[i][j],vectorMag3D(oKO[i][j]));
        }
    }

}

/* ---------------------------------------------------------------------- */

void STLtri::generate_random(double *pos,RanPark *random)
{
    // step 1 - choose triangle
    int chosen = 0;
    double Area_accum = 0;
    double rd = random->uniform() * Area_total;

    while(Area_accum < rd && chosen < nTri)
        Area_accum += Area[chosen++];

    chosen--;

    if(chosen >=nTri) error->all("STLtri::generate_random error");

    // step 2 - random bary coords
    double u = random->uniform();
    double v = random->uniform();

    double tmp = sqrt(u);

    double bary_0 = 1 - tmp;
    double bary_1 = v * tmp;
    double bary_2 = 1 - bary_0 - bary_1;

    // step 3 - apply bary coords
    pos[0] = bary_0 * node[chosen][0][0] + bary_1 * node[chosen][1][0] + bary_2 * node[chosen][2][0];
    pos[1] = bary_0 * node[chosen][0][1] + bary_1 * node[chosen][1][1] + bary_2 * node[chosen][2][1];
    pos[2] = bary_0 * node[chosen][0][2] + bary_1 * node[chosen][1][2] + bary_2 * node[chosen][2][2];

}

/* ---------------------------------------------------------------------- */

void STLtri::getTriBoundBox(int iTri, double *xmin, double *xmax,double delta)
{
    
    xmin[0] = xmin[1] = xmin[2] = +BIG;
    xmax[0] = xmax[1] = xmax[2] = -BIG;

    for(int j=0;j<3;j++)
    {
      xmin[j] = MathExtraLiggghts::min(node[iTri][0][j], MathExtraLiggghts::min(node[iTri][1][j],node[iTri][2][j])) - delta;
      xmax[j] = MathExtraLiggghts::max(node[iTri][0][j], MathExtraLiggghts::max(node[iTri][1][j],node[iTri][2][j])) + delta;
    }
}

/* ---------------------------------------------------------------------- */

void STLtri::getMeshBoundBox(double *xmin, double *xmax)
{
    double xmin_i[3],xmax_i[3];
    
    xmin[0] = xmin[1] = xmin[2] = +BIG;
    xmax[0] = xmax[1] = xmax[2] = -BIG;

    for(int i = 0; i < nTri; i++)
    {
        getTriBoundBox(i, xmin_i, xmax_i, 0.);
        vectorComponentMin3D(xmin,xmin_i,xmin);
        vectorComponentMax3D(xmax,xmax_i,xmax);
    }
}

/* ---------------------------------------------------------------------- */

bool STLtri::is_planar()
{
    bool isPlanar = true;

    for (int iTri = 0; iTri < nTri; iTri++)
    {
        if(!isPlanar) break;
        for(int j = 0; j < 3; j++)
        {
            int jTri = neighfaces[iTri][j];
            if(jTri < 0) continue;
            isPlanar = isPlanar && are_coplanar_neighs(iTri,jTri);
        }
    }

    return isPlanar;
}

/* ---------------------------------------------------------------------- */

bool STLtri::are_coplanar_neighs(int i1,int i2)
{
    
    bool are_neighs;

    if     (neighfaces[i1][0] == i2)
        are_neighs = contactInactive[i1] & EDGE_INACTIVE[0];
    else if(neighfaces[i1][1] == i2)
        are_neighs =  contactInactive[i1] & EDGE_INACTIVE[1];
    else if(neighfaces[i1][2] == i2)
        are_neighs =  contactInactive[i1] & EDGE_INACTIVE[2];
    else return false;

    if     (neighfaces[i2][0] == i1)
        are_neighs = are_neighs && contactInactive[i2] & EDGE_INACTIVE[0];
    else if(neighfaces[i2][1] == i1)
        are_neighs = are_neighs && contactInactive[i2] & EDGE_INACTIVE[1];
    else if(neighfaces[i2][2] == i1)
        are_neighs = are_neighs && contactInactive[i2] & EDGE_INACTIVE[2];
    else return false;

    return are_neighs;
}

/* ---------------------------------------------------------------------- */

void STLtri::normal_vec(double *vec,int i)
{
    vectorCopy3D(facenormal[i],vec);
}

/* ---------------------------------------------------------------------- */

int STLtri::is_on_surface(double *pos)
{
    int on_surf = 0;

    //brute force
    for(int iTri = 0; iTri < nTri; iTri++)
    {
        on_surf += is_in_tri(pos,iTri);
    }

    if(on_surf >= 1) return 1;
    return 0;
}

/* ---------------------------------------------------------------------- */

int STLtri::is_in_tri(double *pos,int i)
{
    double v0[3],v1[3],v2[3];
    double dot00,dot01,dot02,dot11,dot12,invDenom,u,v;

    vectorSubtract3D(node[i][2], node[i][0], v0);
    vectorSubtract3D(node[i][1], node[i][0], v1);
    vectorSubtract3D(pos,        node[i][0], v2);

    dot00 = vectorDot3D(v0, v0);
    dot01 = vectorDot3D(v0, v1);
    dot02 = vectorDot3D(v0, v2);
    dot11 = vectorDot3D(v1, v1);
    dot12 = vectorDot3D(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;

    if((u > -SMALL) && (v > -SMALL) && (u + v < 1+SMALL)) return 1;
    return 0;
}

