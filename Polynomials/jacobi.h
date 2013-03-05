/*
 *  jacobi.h
 *
 *	Fornisce il polinomio di jacobi, la derivata prima, la normalizzazione, le funzioni normalizzate
 *
 *  Jacobi polynomials and Normalization
 *
 *  P_n^(alpha,beta)(x)
 *
 *
 *  Derivatives P'_n^(alpha,beta)(x)
 *
 *  
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

#ifndef JAC_POLI_H
#define JAC_POLI_H

#include "algebra.h"
#include <assert.h>


//Polinomio di jacobi P_n^(alpha,beta)(x)
inline long double jacobi(const unsigned short & n, 
						  const long double & alpha, 
						  const long double & beta, 
						  const long double & x) {
	assert (alpha >= -1. && beta >= -1. && std::fabs(x) <= 1. );						
	if ( n == 0 ) return 1;																					
	else if ( n == 1 ) return  ( ( alpha - beta + (alpha + beta + 2. ) * x ) / 2. ) ;								

	else {
		long double P[ n + 1 ];
		P[0] = 1.;
		P[1] = ( alpha - beta + (alpha + beta + 2. ) * x ) / 2.;
		for (short i = 1 ; i < n  ; ++i) {
			P [ i + 1 ] = ( 
							( (2. * i + alpha + beta + 1.) * ( pow( alpha , 2 ) - pow( beta , 2 ) ) + pochhammer( (long double) 2. * i + alpha + beta , 3. ) * x ) * P[i]
							- 2. * ( i + alpha ) * ( i + beta ) * ( 2. * i + alpha + beta + 2 ) * P[i-1]
						   )
							/ 
						( 2. * (i + 1.) * ( i + alpha + beta + 1.) * (2*i + alpha + beta) ) ;
		}
		return P[ n ];
	}
	
	
	
// 	1. - jacobi (n1 , alpha , beta ,  1. )/c_binomiale ( n1 + alpha , n1) == 0
}

//P'_n^(alpha,beta)(x)
inline long double djacobi(const unsigned short & n, 
						   const long double & alpha, 
						   const long double & beta, 
						   const long double & x) {
		
	assert (alpha >= -1. && beta >= -1. && std::fabs(x) <= 1. );						
	if ( n == 0 ) return 0;																					
	else if ( n == 1 ) return  ( (alpha + beta + 2. ) / 2. ) ;														
	else return ( ( n + alpha + beta + 1. ) * jacobi ( n - 1 , alpha + 1. , beta + 1. , x ) / 2. );
}

//normalizzazione per il polinomo di jacobi
inline long double N_jacobi(const unsigned short & n, 
							const long double & alpha, 
							const long double & beta) {

	return sqrt(
				(2*n + alpha + beta + 1)					* pochhammer ( n + beta + 1., alpha)  / 
				( pow((long double)2., alpha + beta + 1)	* pochhammer ( n + 1., alpha)  )
				);
	
}

	
//Normalized P_n^(alpha,beta)(x)
inline long double jacobi_N(const unsigned short & n,
							const long double & alpha, 
							const long double & beta, 
							const long double & x) {
	return  N_jacobi(n,alpha,beta) * jacobi(n,alpha,beta,x);
}

//Normalized P'_n^(alpha,beta)(x)
inline long double djacobi_N(const unsigned short & n,
							 const long double & alpha, 
							 const long double & beta, 
							 const long double & x) {
	return  N_jacobi(n,alpha,beta) * djacobi(n,alpha,beta,x);
}


#endif
