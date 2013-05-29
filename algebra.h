/*
 *  algebra.h
 *  Tesi
 *
 *  Created by Marco Resa on 23/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef ALGEBRA_H
#define ALGEBRA_H

#include <cmath>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/beta.hpp>
#include "costanti_fisiche.h"


//definisce il tipo complesso


//Delta di Kronecker
inline bool Kdelta (const int & n, const int & m ) {
	return ( n == m );
}

// elevamento di - 1 a potenza intera
inline short menouno (const int & n ) {
	return (n %2) ? -1 : 1;
//	if ( n %2 ) return -1; // se pari n%2 = 0 = falso
//	else return 1;
}

inline short segno (const long double & n ) {
	return (n < 0.) ? -1 : 1;
//	if ( n %2 ) return -1; // se pari n%2 = 0 = falso
//	else return 1;
}

//simbolo di pochhammer = (a)_b = tgamma(a+b)/tgamma(a)
inline long double pochhammer (const long double &a , 
							   const long double &b ) {
	return 1./boost::math::tgamma_delta_ratio(a, b);
}

//coefficiente binomiale
inline long double c_binomiale (const long double &a , 
								const long double &b ) {
	return 1./( b * boost::math::beta(b, a - b + 1.));
}

#endif /* ALGEBRA_H */
