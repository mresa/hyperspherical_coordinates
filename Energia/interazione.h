/*
 *  interazione.h
 *  Tesi
 *
 *  Created by Marco Resa on 17/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef INTERAZIONE_H
#define INTERAZIONE_H
#include <cmath>

inline long double interazione (const long double & r){
	/*long double E_1 = 144.86 ;
	long double R_1 = 0.82 ;
	long double E_2 = -83.34 ;
	long double R_2 = 1.6 ;
	return E_1 * exp ( - pow( r / R_1 , 2) ) + E_2 * exp ( - pow( r / R_2 , 2) );
	 */
	long double E_1 = 144.86 ;
	long double R_1 = 0.82 ;
	long double E_2 = -83.34 ;
	long double R_2 = 1.6 ;
	return E_1 * exp ( - pow( r / R_1 , 2) ) + E_2 * exp ( - pow( r / R_2 , 2) );
	//return 1.;
}

inline long double interazione_coulombiana (const long double & r){
	//return struttura_fine * h_bar_c / r ;
	return e_quadro / r ;
}


//spinoriale

inline long double interazione_0 (const long double & r){
	long double E_1 = 144.86 ;
	long double R_1 = 0.82 ;
	long double E_2 = -66.7 ;
	long double R_2 = 1.6 ;
	//
	return E_1 * exp ( - pow( r / R_1 , 2) ) + E_2 * exp ( - pow( r / R_2 , 2) );
	
	//return interazione(r);
}

inline long double interazione_1 (const long double & r){
	long double E_1 = 144.86 ;
	long double R_1 = 0.82 ;
	long double E_2 = -97.0 ;
	long double R_2 = 1.6 ;
	//
	return E_1 * exp ( - pow( r / R_1 , 2) ) + E_2 * exp ( - pow( r / R_2 , 2) );
	
	//return interazione(r);

}


inline long double interazione_s_0 (const long double & r){
	return (interazione_0(r) + 3.*interazione_1(r))/4.;
}

inline long double interazione_s_1 (const long double & r){
	return (-interazione_0(r) + interazione_1(r))/4.;
}

//isospinoriale
/* sbagliata
inline long double interazione_st_00 (const long double & r){
	return (interazione_0(r) + 9.*interazione_1(r))/16.;
	//	return interazione(r);
}

inline long double interazione_st_01 (const long double & r){
	return (-interazione_0(r) + 3.*interazione_1(r))/16.;
	//	return 0.;
}

inline long double interazione_st_10 (const long double & r){
	return (-interazione_0(r) + 3.*interazione_1(r))/16.;
	//	return 0.;
}

inline long double interazione_st_11 (const long double & r){
	return (interazione_0(r) + interazione_1(r))/16.;
	//	return 0.;
}
*/
//giusta
inline long double interazione_st_00 (const long double & r){
	return (interazione_0(r) + interazione_1(r))*3./16.;
//	return interazione(r);
}

inline long double interazione_st_01 (const long double & r){
	return (interazione_0(r) -3.*interazione_1(r))/16.;
//	return 0.;
}

inline long double interazione_st_10 (const long double & r){
	return (interazione_1(r) - 3.*interazione_0(r))/16.;
//	return 0.;
}

inline long double interazione_st_11 (const long double & r){
	return -(interazione_0(r) + interazione_1(r))/16.;
//	return 0.;
}
#endif
