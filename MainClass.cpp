/*
 * MainClass.cpp
 *
 *  Created on: Mar 1, 2012
 *      Author: jn
 */
#include <iostream>
#include <math.h>

#define PI 		atan(1) * 4
#define NB_INT 	100
#define RADIUS 	5
#define DIST 	10 * RADIUS

using namespace std;

struct param {
  float da;
  float db;
  float D;

  float cosa;
  float cosb;
  float sina;
  float sinb;

  float R [2];
  float integral;

};

//To extend to a third D integral. (integral depending only on one dimension)
float third_D_coef(float d, float sind)
{
	return 2 * PI * d * sind;
}

//returns the squared norm of the vector vect (2 Dimensions)
float sqr_vect_norm(float * vect)
{
	return pow(vect[0], 2) + pow(vect[1], 2);
}

//Calculates the contribution of one element in the integral
float integrate(struct param * par)
{
	return (par->cosa * par->cosb / sqr_vect_norm(par->R));
}

void calculate_vector(struct param * par, int ia, int ib)
{
	//Calculates the value of the R vector
	par->R[0] = par->D - par->cosa - par->cosb;
	par->R[1] = - par->sina - par->sinb;
}

//This function tests if the two points can physically "see" each other,
//by testing the scalar product of the 2 vectors.
bool scalar_test(struct param * par)
{
	bool test = false;
	float scala = par->cosa * par->R[0] + par->sina * par->R[1];
	float scalb = (-par->cosb) * par->R[0] + (-par->sinb) * par->R[1];

	if(scala > 0)
		test = true;

	if(scalb > 0)
		test = true;

	return test;
}

int procedeCalc(float rad, float dist)
{
	int ia = 0;
	int ib = 0;
	float integral = 0;

	//init parameter

	struct param par;

	//If the radius is not specified, take the default value
	if ( rad >= 0)
		par.da = par.db = rad;
	else
		par.da = par.db = RADIUS;

	//If the distance is not specified, take the default value
	if( dist >= 0)
		par.D = dist;
	else
		par.D = DIST;

	par.integral = 0;

	int nb = NB_INT; //Number if integration steps

	//starting discrete integration
	for (ia=0;ia<nb;ia++)
	{
		integral = 0;

		// Calculates the vector da coordinates, according to the angle ia
		par.cosa = par.da * cos((PI * ia) / (2 * NB_INT));
		par.sina = par.da * sin((PI * ia) / (2 * NB_INT));

		for (ib=0;ib<nb;ib++)
		{

			// Calculates the vector db coordinates, according to the angle ib
			par.cosb = par.db * cos((PI * ib) / (2 * NB_INT));
			par.sinb = par.db * sin((PI * ib) / (2 * NB_INT));

			calculate_vector(&par, ia, ib);
			if(scalar_test(&par))
				integral = integral + integrate(&par) * third_D_coef(par.db, par.sinb);
		}
		integral = integral * third_D_coef(par.da, par.sina);
		par.integral = par.integral + integral;
	}

	par.integral = integral / (PI * 2 * PI * pow(par.da, 2)); //normalize integral
	cout 	<< "distance : " << par.D
			<< ", radius : " << par.da
			<< ", shape factor = " << par.integral
			<< endl;

	return 0;
}

int main()
{
	int r = 0;
	cout << "Hello World " << PI << endl ;
	for (int rad = 1; rad < 5 ; rad++)
	{
		for (int dist = 2 * rad; dist<50; dist++)
		{
			//provide the distance and radius value, or -1
			//if -1, default values are taken.
			//default radius = 5
			//default distance = 10 * RADIUS
			r = procedeCalc(rad, dist);
		}
	}
	return 0;
}

