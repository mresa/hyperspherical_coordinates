/*
 *  integrale_angolare.h
 *  Tesi
 *
 *  Created by Marco Resa on 21/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef INTEGRALE_ANGOLARE_H
#define INTEGRALE_ANGOLARE_H
#include "clebsch.h"

//integrale di 3 armoniche sferiche, ridotto  -> < j1,0;j2,0|j,0>
inline long double int_3j_ridotto (const unsigned short& j1, 
								   const unsigned short& j2, 
								   const unsigned short& j) {
	if (dis_triang(j, j1, j2)) return sqrt((long double)(2.*j1+1)*(2.*j2+1)/(4*pi*(2*j + 1))) * clebsch_gordan(j1,j2,j);
	else return 0.;
}

//integrale di 3 armoniche sferiche, parte in m
inline long double int_3j_m (const unsigned short& j1, 
							 const unsigned short& j2, 
							 const unsigned short& j,
							 const short& m1, 
							 const short& m2, 
							 const short& m) {
	if (dis_triang(j, j1, j2) && m1 + m2 == m) return clebsch_gordan(j1,j2,j,m1,m2,m);
	else return 0.;
}

//integrale di 3 armoniche sferiche
inline long double int_3j (const unsigned short& j1, 
						   const unsigned short& j2, 
						   const unsigned short& j,
						   const short& m1, 
						   const short& m2, 
						   const short& m) {
	long double risultato = int_3j_m(j1,j2,j,m1,m2,m);
	if (risultato !=0) return int_3j_ridotto(j1,j2,j)*risultato;
	else return 0.;
}

//integrale di 2 armoniche sferiche
inline long double int_2j (const unsigned short& j1, 
						   const unsigned short& j2, 
						   const short& m1, 
						   const short& m2) {
	return sqrt((long double) 4*pi) * int_3j(j1,0,j2,m1,0,m2);
}

#endif

