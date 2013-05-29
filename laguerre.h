/*
 *  laguerre.h
 *  
 *  Laguerre polynomials
 *  
 *  L_n(x) == laguerre(n,x)
 *  L_n^k(x) == laguerre(n,k,x)
 *  
 *  and derivatives
 *  
 *  L'_n(x) == dlaguerre(n,x)
 *  L'_n^k(x) == dlaguerre(n,k,x)
 *  
 */

/*
 *	Copyright (C) 2013 by Marco Resa <https://github.com/mresa>
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef LAG_POLI_H
#define LAG_POLI_H

#include <cmath>
#include <gsl/gsl_sf_bessel.h>
#include <boost/math/special_functions/laguerre.hpp>
#include <assert.h>

//Polinomio di laguerre L_n^alpha(x)
inline long double laguerre(const unsigned short & n, 
							const short & alpha, 
							const long double & x) {
	
	assert(alpha > -1 && x >= 0 );												// controllo dei limiti: alpha  > -1
	if (n==0) return 1;															// L_0^alpha(x) == 1
	else if (n==1) return (- x + alpha + 1.);									// L_1^alpha(x) == alpha + 1 - x
	else return boost::math::laguerre(n, alpha , x);							// utilizzo la boostlibrary
	
}

//L_n(x)
inline long double laguerre(const unsigned short & n,
							const long double & x) {
	return laguerre(n, 0, x);
}

//L'_n^alpha(x)
inline long double dlaguerre(const unsigned short & n, 
							 const short & alpha, 
							 const long double & x) {
	
	assert(alpha > -1 && x >= 0 );                                                      // controllo dei limiti: alpha  > -1
	if (n==0) return 0;                                                                 // L'_0 == 0
	else if (n==1) return (-1.);														// L'_1(x) == -1
	else return (n*laguerre(n, alpha, x) - (n + alpha)*laguerre(n-1, alpha, x) ) /x;	// calcolo iterativo di L'_l(x)
}

//L'_n(x)
inline long double dlaguerre(const unsigned short & n,
							 const long double & x) {
	return dlaguerre(n, 0, x);
}

//funzione di controllo per gli zeri ottenuti
//check zeros
inline bool controllo_zeri_laguerre(const unsigned short & n, 
									const short & alpha, 
									const long double & x, 
									const unsigned short & m) {
	
	//limiti per lo zero m-esimo di L_n, secondo abramovitz-stegun	
	assert(alpha > -1 && x >= 0 );												// controllo dei limiti: alpha  > -1
	float k_n = n + (1. + alpha)/2.;	
	float k_m = m + (1. + alpha)/2.;	
	return ( x > pow(gsl_sf_bessel_zero_Jnu(alpha, m),2)/(4*k_n) && x < (k_m/k_n)*( 2.*k_m + sqrt(4*(k_m*k_m)+ 1/4) ) );
}

inline bool controllo_zeri_laguerre(const unsigned short & n, 
									const long double & x, 
									const unsigned short & m) {
	return controllo_zeri_laguerre(n, 0, x, m);
}

//prima approssimazione per gli zeri, zero m-esimo
inline long double guess_zero(const unsigned short & n, 
							  const short & alpha, 
							  const unsigned short & m) {
/*
 *  Guess iniziale per gli zeri: referenze "Numerical recipes in fotran77 ...", Cambridge U Press, 1999 oppure
 *  <http://people.sc.fsu.edu/~jburkardt/cpp_src/laguerre_rule/laguerre_rule.cpp>
 */
	assert(alpha > -1);												// controllo dei limiti: alpha  > -1
	long double zeri[m];
	for (short i = 0 ; i < m ; ++i) { 
		if ( i == 0 ) {
			zeri[i] = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
			( 1.0 + 2.4 * n + 1.8 * alpha );
		}
		else if ( i == 1 ) {
			zeri[i] = zeri[i-1] + ( 15.0 + 6.25 * alpha ) / 
			( 1.0 + 0.9 * alpha + 2.5 * n );
		}
		else {
			long double r1 = ( 1.0 + 2.55 * ( long double ) ( i - 1 ) ) 
				/ ( 1.9 * (long double ) ( i - 1 ) );
				
			long double r2 = 1.26 * ( long double ) ( i - 1 ) * alpha / 
				( 1.0 + 3.5 * ( long double ) ( i - 1 ) );
				
			long double ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );
				
			zeri[i] = zeri[i-1] + ratio * ( zeri[i-1] - zeri[i-2] );
		}
	}
	for(short i = 0 ; i<1000;++i ) {
		zeri[m-1] -= (laguerre(n, alpha, zeri[m-1])/dlaguerre(n, alpha, zeri[m-1]));
	}
	return zeri[m-1];
}

#endif
