/*
 *  gauss-legendre.h
 *  Tesi
 *
 *  Created by Marco Resa on 23/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef GAUSS_LEG_H
#define GAUSS_LEG_H
#include "err.h"
#include "legendre.h"
#include "algebra.h"
#include <assert.h>
#include <omp.h> //calcolo parallelo OMP

/*
 *  
 *  
 *	Punti  e pesi per la formula 25.4.29 integrazione tra 1 e -1 
 *
 *
 */

class Gauss_Legendre {
	unsigned short n;						//ordine laguerre
	long double* x;							//griglia degli zeri
	long double* w;							//griglia dei pesi
	
	static inline bool errore_gauss(){
		std::cerr << "errore in gauss-legendre.h" << std::endl;
		std::exit(1);
		return 1;
	}
	
	static inline bool errore_zeri(){
		std::cerr << "errore negli zeri" << std::endl;
		std::exit(1);
		return 1;
	}
	
	//controllo per gli zeri dei polinomi di legendre
	static inline bool controllo(const unsigned short & n, 
								 const unsigned short & i, 
								 const long double & x) { //controllo dei limiti per lo zero, per vedere se ho preso quello giusto
		return ( std::acos(x) >= (2*i+1)*pi/(2*n+1) && std::acos(x) <= 2*(i+1)*pi/(2*n+1) );
	}
	
	static inline void nodi(const unsigned short & n, 
							long double* zeri, 
							long double* pesi) {
		
		if ( n<1 ) {													// controllo dei limiti: L_n(x) ha n zeri; in particolare L_0 non ne ha
			errore_gauss();
		}
		else {
			
			/*
			 *	 Calcolo degli zeri di legendre: guess iniziale e algoritmo di newton
			 */
			
			
		#pragma omp parallel for schedule( guided ) //ciclo for parallelo, senza membri comuni o privati
			for (int i = 0 ; i < n ; ++i) { 
				zeri[i] = cos( (4*i + 3 ) * pi / (4*n+2) + 1/(8*n*n*std::tan( (4*i + 3) * pi / (4*n-2) ) ) );
				
				/*
				 *	Algoritmo newton per il calcolo degli zeri: si ferma per incrementi minori di err
				 */
				
				while (std::fabs(legendre(n, zeri[i])/dlegendre(n, zeri[i])) > (err*zeri[i]) && std::fabs(legendre(n, zeri[i])/dlegendre(n, zeri[i])) > err){
					zeri[i] -= (legendre(n, zeri[i])/dlegendre(n, zeri[i]));
				}
				
				
				/* 
				 *  controllo dei limiti per lo zero trovato
				 */
				
				//				if ( controllo(n, i , zeri[i]));
				//				else errore_zeri() ;
				
				pesi[i]=2./((1 - pow( (zeri[i]) , 2 )) * pow( dlegendre( n, zeri[i] ), 2 ));
			}
		}
	}
	
public:
	
	
	Gauss_Legendre ( const unsigned short & N ) {
		if (N <= 0) errore_gauss();
		else {
			n = N;
			x = new long double[n];
			w = new long double[n];
			nodi(n, x, w);
		}	
	}	
	
	~Gauss_Legendre () {
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

/*
 *  
 *  
 *	Punti  e pesi per la formula 25.4.35 integrazione tra 1 e -1 
 *
 *
 */

class Gauss_Legendre_sqrt {
	unsigned short n;						//ordine laguerre
	long double* x;							//griglia degli zeri
	long double* w;							//griglia dei pesi
	
	static bool errore_gauss(){
		std::cerr << "errore in gauss-legendre_sqrt.h" << std::endl;
		std::exit(1);
		return 1;
	}
		
public:
	
	
	Gauss_Legendre_sqrt (const unsigned short N) {
		if (N <= 0) errore_gauss();
		else {
			n = N;
			x = new long double[n];
			w = new long double[n];
			unsigned short n_G = 2*n+1;
			Gauss_Legendre G(n_G);
			unsigned short n_controllo = 0;
			for (unsigned short i = 0; i < n_G && n_controllo < n ; ++i) {
				if (G.X(i) > 0.) {
					x[n_controllo] = -1. + 2. * (1.-pow( G.X(i), 2));
					w[n_controllo] = 4. * sqrt( 2. ) * pow( G.X(i), 2) * G.W(i);
					++n_controllo;
				}
			}
			if (n_controllo != n) {
				errore_gauss();
			}
		}	
	}	
	
	
	
	~Gauss_Legendre_sqrt () {
		delete[] x;
		delete[] w;
	}
	
	long double X(const unsigned short & i) const {
		assert( i < n );
		return x[i];
	}
	long double W(const unsigned short & i) const {
		assert( i < n );
		return w[i];
	}
	unsigned short N() const {return n;}
};

#endif
