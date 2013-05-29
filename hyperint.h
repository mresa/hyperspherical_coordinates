/*
 *  HyperInt.h
 *  Tesi
 *
 *  Created by Marco Resa on 28/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
/*	pesi e punti per l'integrazione dz (1-z)^dubalpha/2 (1+z)^dubbeta/2 f(z)
 *
 *	con dubalpha e dubbeta interi
 */

#ifndef HYPERINT_H
#define HYPERINT_H

#include <assert.h>
#include "hyper.h"
#include "algebra.h"	
#include "clebsch.h"	
#include "gauss-chebyshev.h"
#include "gauss-legendre.h"
#include <omp.h> //calcolo parallelo OMP


class HyperInt {
	unsigned short n;						//ordine integrazione
	long double* x_pp;						//griglia degli zeri
	long double* w_pp;						//griglia dei pesi
	long double* x_dp;						//griglia degli zeri
	long double* w_dp;						//griglia dei pesi
	long double* x_pd;						//griglia degli zeri
	long double* w_pd;						//griglia dei pesi
	long double* x_dd;						//griglia degli zeri
	long double* w_dd;						//griglia dei pesi
	
public:
	
	
	HyperInt (const unsigned short & N) {
		assert (N > 0);
		n = N;
		x_pp = new long double[n];
		w_pp = new long double[n];
		x_dp = new long double[n];
		w_dp = new long double[n];
		x_pd = new long double[n];
		w_pd = new long double[n];
		x_dd = new long double[n];
		w_dd = new long double[n];
		
		Gauss_Legendre G (n);			
#pragma omp parallel for schedule( guided )
		for (short i = 0 ; i < n; ++i) {
			x_pp[i] = G.X(i);
			w_pp[i] = G.W(i);
		}
		
		Gauss_Legendre_sqrt G1 (n);			
#pragma omp parallel for schedule( guided )
		for (short i = 0 ; i < n; ++i) {
			x_pd[i] = -G1.X(i);
			w_pd[i] = G1.W(i);
	
			x_dp[i] = G1.X(i);
			w_dp[i] = G1.W(i);
		}						
		
		Gauss_Chebyshev G2 (n);			
#pragma omp parallel for schedule( guided )
		for (short i = 0 ; i < n; ++i) {
			x_dd[i] = G2.X(i);
			w_dd[i] = G2.W(i);
		}	
	}
	
	~HyperInt () {
		delete[] x_pp;
		delete[] w_pp;
		delete[] x_pd;
		delete[] w_pd;
		delete[] x_dp;
		delete[] w_dp;
		delete[] x_dd;
		delete[] w_dd;
	}
	
	inline long double X(const unsigned short & i, 
						 const unsigned short & DUBALPHA, 
						 const unsigned short & DUBBETA) const {
		assert( i < n );
		if(DUBALPHA %2 == 0){
			if(DUBBETA %2 == 0){
				return x_pp[i];
			}
			else{ //legendre 2 modificata
				return x_pd[i];
			}						
		}
		else {
			// alpha dispari, alpha/2 semintero
			if(DUBBETA %2 == 0){  // legendre 2 non modificata
				return x_dp[i];
			}
			
			// alpha semintero, beta semintero
			else{
				return x_dd[i];
			}
		}
	}
	
	inline long double W(const unsigned short & i, 
						 const unsigned short & DUBALPHA, 
						 const unsigned short & DUBBETA) const {
		assert( i < n );
		if(DUBALPHA %2 == 0){
			if(DUBBETA %2 == 0){
				return w_pp[i];
			}
			else{ //legendre modificata
				return w_pd[i];
			}						
		}
		else {
			// alpha dispari, alpha/2 semintero
			if(DUBBETA %2 == 0){  // legendre 2 non modificata
				return w_dp[i];
			}
			
			// alpha semintero, beta semintero
			else{
				return w_dd[i];
			}
		}
	}
	
	inline long double N() const {return n;}
};

/*	pesi e punti per l'integrazione dz (1-z)^dubalpha/2 (1+z)^dubbeta/2 f(z)
 *
 *	con dubalpha e dubbeta interi
 *
 * con questa classe basta fare sum( f(x_i) * w_i : pensa lui al resto
 */
class Integrale_Pippo {
	unsigned short n;						//ordine integrazione
	unsigned short alpha;						
	unsigned short beta;						
	long double* x;							//griglia degli zeri
	long double* w;							//griglia dei pesi
	
public:
	
	
	Integrale_Pippo (const HyperInt & G, 
					 const unsigned short & DUBALPHA, 
					 const unsigned short & DUBBETA) {
		n = G.N();
		assert (n > 0);
		alpha = DUBALPHA;
		beta = DUBBETA;
		x = new long double[n];
		w = new long double[n];
		
		if(alpha %2 == 0){
			if(beta %2 == 0){
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i,alpha,beta);
					w[i] = (pow(1. - G.X(i,alpha,beta) , (long double) alpha/2.)	*
							pow(1. + G.X(i,alpha,beta) , (long double) beta/2.)	* 
							G.W(i,alpha,beta) );
				}
			}
			else{ //legendre modificata
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i,alpha,beta);
					w[i] = (pow(1. - G.X(i,alpha,beta) , (long double) alpha/2.)	*
							pow(1. + G.X(i,alpha,beta) , (long double) beta/2. - 1./2.)	* 
							G.W(i,alpha,beta) );
					
				}						
			}
		}
		
		else {
			// alpha dispari, alpha/2 semintero
			if(beta %2 == 0){  // legendre 2 non modificata
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i,alpha,beta);
					w[i] = (pow(1. - G.X(i,alpha,beta) , (long double) alpha/2. - 1./2.)	*
							pow(1. + G.X(i,alpha,beta) , (long double) beta/2.)	* 
							G.W(i,alpha,beta) );
					
				}
			}
			
			// alpha semintero, beta semintero
			else{
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i,alpha,beta);
					w[i] = (pow( 1. - G.X(i,alpha,beta) , (long double) alpha/2. - 1./2.)	*
							pow(1. + G.X(i,alpha,beta) , (long double) beta/2. - 1./2.)	* 
							G.W(i,alpha,beta) );
				}						
			}
		}
	}
	
/*	Integrale_Pippo (const unsigned short & N, 
					 const unsigned short & DUBALPHA, 
					 const unsigned short & DUBBETA ) {
		assert (N > 0);
		n = N;
		alpha = DUBALPHA;
		beta = DUBBETA;
		x = new long double[n];
		w = new long double[n];
		
		if(alpha %2 == 0){
			if(beta %2 == 0){
				Gauss_Legendre G(n);			
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i);
					w[i] = (pow(1. - G.X(i) , (long double) alpha/2.)	*
							pow(1. + G.X(i) , (long double) beta/2.)	* 
							G.W(i) );
				}
			}
			else{ //legendre modificata
				Gauss_Legendre_sqrt G(n);			
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = -G.X(i);
					w[i] = (pow(1. + G.X(i) , (long double) alpha/2.)	*
							pow(1. - G.X(i) , (long double) (beta - 1)/2.)	* 
							G.W(i) );
					
				}						
			}
		}
		
		else {
			// alpha dispari, alpha/2 semintero
			if(beta %2 == 0){  // legendre 2 non modificata
				Gauss_Legendre_sqrt G(n);			
#pragma omp parallel for schedule( guided )
				for (short i = 0 ; i < n; ++i) {
					x[i] = G.X(i);
					w[i] = (pow(1. - G.X(i) , (long double) (alpha - 1)/2.)	*
							pow(1. + G.X(i) , (long double) beta/2.)	* 
							G.W(i) );
					
				}
			}
			
			// alpha semintero, beta semintero
			else{
				Gauss_Chebyshev G(n);							
				if (alpha > beta){
#pragma omp parallel for schedule( guided )
					for (short i = 0 ; i < n; ++i) {
						x[i] = G.X(i);
						w[i] = (pow(1. - G.X(i) * G.X(i), (long double) (beta - 1)/2.)	*
								pow(1. - G.X(i) , (long double) (alpha - beta)/2.)	* 
								G.W(i) );
					}						
				}
				else{
#pragma omp parallel for schedule( guided )
					for (short i = 0 ; i < n; ++i) {
						x[i] = G.X(i);
						w[i] = (pow(1. - G.X(i) * G.X(i), (long double) (alpha - 1)/2.)	*
								pow(1. + G.X(i) , (long double) (beta - alpha )/2.)	* 
								G.W(i) );
					}						
				}
				
				
			}
		}
	}
*/	
	~Integrale_Pippo () {
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
	inline long double ALPHA() const {return alpha;}
	inline long double BETA() const {return beta;}
};


#endif

