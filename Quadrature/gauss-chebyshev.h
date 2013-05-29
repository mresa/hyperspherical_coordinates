/*
 *  gauss-chebyshev.h
 *
 *  Tesi
 *
 *	Punti  e pesi per la formula 25.4.40 integrazione tra 1 e -1, peso sqrt(1-x^2)
 *
 *
 *  Created by Marco Resa on 04/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GAUSS_CHEB_H
#define GAUSS_CHEB_H
#include "err.h"
#include "algebra.h"
#include <assert.h>
#include <omp.h> //calcolo parallelo OMP

class Gauss_Chebyshev {
	unsigned short n;						//ordine integrazione
	long double* x;							//griglia degli zeri
	long double* w;							//griglia dei pesi
	
	static inline bool errore_gauss(){
//		std::cerr << "errore in gauss-chebyshev.h" << std::endl;
//		std::exit(1);
		return 1;
	}
		
	static inline void init(const unsigned short & n, 
							long double* zeri, 
							long double* pesi) {
		
		#pragma omp parallel for schedule( guided ) //ciclo for parallelo, senza membri comuni o privati
		for ( short i = 0 ; i < n ; ++i) { 
			long double I = ( i + 1. ) / ( n + 1. ) * pi ;
			zeri[i] = std::cos( I );
			pesi[i] = pi * pow( std::sin( I ) , 2 ) / ( n + 1. );
		}
	}
	
public:
	
	
	Gauss_Chebyshev ( const unsigned short & N ) {
		if (N < 1) errore_gauss();
		else {
			n = N;
			x = new long double[n];
			w = new long double[n];
			init(n , x , w);
		}	
	}
	
	~Gauss_Chebyshev () {
		delete[] x;
		delete[] w;
	}
	
	inline long double X(const unsigned short & i) const {
		assert( i < n );
		return x[i];
	}
	inline long double W(const unsigned short & i) const {
		assert( i < n );
		return w[i];
	}
	inline unsigned short N() const {return n;}
};


#endif
