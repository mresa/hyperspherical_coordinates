/*
 *  clebsh.h
 *  Tesi
 *	Fornisce :
 	
 *	i coefficienti di clebsh gordan per < j1,m1;j2,m2|j,m> , in doppia precisione
 
 *	il risultato di una integrazione di 2 o 3 armoniche sferiche

 *  Created by Marco Resa on 13/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CLEBSCH_H
#define CLEBSCH_H
#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "algebra.h"
#include "algebra_vettori.h"
#include "triangolari.h"

// < j1,m1;j2,m2|j,m>
inline long double clebsch_gordan (const long double& j1, 
								   const long double& j2, 
								   const long double& j,
								   const long double& m1, 
								   const long double& m2, 
								   const long double& m) {
 

	return (long double) (
						  menouno(m+j1-j2)
						  *
						  sqrt( (long double) 2.*j + 1. )
						  *
						  gsl_sf_coupling_3j((int) ( 2*j1 ), 
											 (int) ( 2*j2 ), 
											 (int) ( 2*j  ),
											 (int) ( 2*m1 ),
											 (int) ( 2*m2 ),
											 (int) ( -2*m )));
}

// < j1,0;j2,0|j,0>
inline long double clebsch_gordan (const long double& j1, 
								   const long double& j2, 
								   const long double& j) {
	return clebsch_gordan (j1, j2, j, 0., 0., 0.);
}


// < j1,m1;j2,m2|j,m>
inline long double clebsch_gordan (const unsigned short& j1, 
								   const unsigned short& j2, 
								   const unsigned short& j,
								   const short& m1, 
								   const short& m2, 
								   const short& m) {
	return (long double) (
						  menouno(m+j1-j2)
						  *
						  sqrt((long double)  2.*j + 1. )
						  *
						  gsl_sf_coupling_3j((int) ( 2*j1 ), 
											 (int) ( 2*j2 ), 
											 (int) ( 2*j  ),
											 (int) ( 2*m1 ),
											 (int) ( 2*m2 ),
											 (int) ( -2*m )));
}

// < j1,0;j2,0|j,0>
inline long double clebsch_gordan (const unsigned short& j1, 
								   const unsigned short& j2, 
								   const unsigned short& j) {
	return clebsch_gordan (j1, j2, j, 0, 0, 0);
}

	
	
#endif

