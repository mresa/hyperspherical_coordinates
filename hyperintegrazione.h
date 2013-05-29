/*
 *  hyperintegrazione.h
 *  Tesi
 *
 *  Created by Marco Resa on 29/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef HYPERINTEGRAZIONE_H
#define HYPERINTEGRAZIONE_H

#include "hyperint.h"

inline bool hyper_Kdelta(const hyper& H1, 
						 const hyper& H2){
	if (H1.N() != H2.N()) return false;
	bool temp = true;
	for (short i = 0; i < H1.N() ; ++i){
		temp *= Kdelta(H1.n(i), H2.n(i));
		temp *= Kdelta(H1.l(i), H2.l(i));
		temp *= Kdelta(H1.m(i), H2.m(i));
	}
	return temp;
}

inline bool hyper_Kdelta_to(const hyper& H1, 
							const hyper& H2,
							unsigned short & indice){
	if (H1.N() != H2.N() || indice >= H1.N()) return false;
	bool temp = true;
	for (short i = 0; i < indice ; ++i){
		temp *= Kdelta(H1.n(i), H2.n(i));
		temp *= Kdelta(H1.l(i), H2.l(i));
		temp *= Kdelta(H1.m(i), H2.m(i));
	}
	return temp;
}


inline bool dis_triang_H_k (const hyper& H2, const hyper& H1, const hyper& H3) {
	bool verita = dis_triang(H2.k(),H1.k(),H3.k());
	return verita;
}

inline bool dis_triang_H_l (const hyper& H2, const hyper& H1, const hyper& H3) {
	bool verita = dis_triang(H2.k(),H1.k(),H3.k()) && dis_triang(H2.l(),H1.l(),H3.l());
	return verita;
}
inline bool dis_triang_H_m (const hyper& H2, const hyper& H1, const hyper& H3) {
	bool verita = dis_triang_m(H2.m(),H1.m(),H3.m());
	return verita;
}
inline bool dis_triang_H (const hyper& H2, const hyper& H1, const hyper& H3) {
	bool verita = dis_triang(H2.k(),H1.k(),H3.k()) && dis_triang(H2.l(),H1.l(),H3.l()) && dis_triang_m(H2.m(),H1.m(),H3.m());
	return verita;
}

//integrale di 2 funzioni di ripelle, all'ordine "ordine", per indice i
inline long double hyper_int_ripelle (const HyperInt & G,
									  const unsigned short& i,
									  const hyper& H1, 
									  const hyper& H2){
	unsigned short ordine = G.N();
	assert (ordine != 0 && i < H1.N() && i < H2.N());
	
	if (i == 0) return 1.;
	
	//controllo se i numeri quantici nu al passo precedente sono uguali
	
	else if (H1.alpha(i) != H2.alpha(i) || H1.beta(i) != H2.beta(i)) return 0.;
	
	//altrimenti se l'indice di particella è dispari (i parte da 0)
	else if ( (i+1) %2 != 0 ) {	//dispari		
		long double integrale = 0.;
		unsigned short n1 = H1.n(i);
		unsigned short n2 = H2.n(i);
		long double alpha = H1.alpha(i);
		long double beta = H1.beta(i);
		
		
		// correzione:  x -> -x
#pragma omp parallel for reduction(+:integrale) schedule( guided )
		for (short i = 0 ; i < ordine; ++i) {
			integrale +=
			jacobi_N( n1, alpha, beta, G.X(i, 2, 1))			*
			jacobi_N( n2, alpha, beta, G.X(i, 2, 1))			*
			pow(1. - G.X(i, 2, 1) , alpha)						* 
			pow(1. + G.X(i, 2, 1) , beta - 1./2.)					* 
			G.W(i, 2, 1);
		}
		integrale /= ( pow((long double)2., alpha + beta + 2));
		return integrale;
	}
	
	else {	//pari
		long double integrale = 0.;
		unsigned short n1 = H1.n(i);
		unsigned short n2 = H2.n(i);
		long double alpha =  H1.alpha(i);
		long double beta = H1.beta(i);
#pragma omp parallel for reduction(+:integrale) schedule( guided )
		for (short i = 0 ; i < ordine; ++i) {
			integrale +=
			jacobi_N( n1, alpha, beta, G.X(i,1,1))			*
			jacobi_N( n2, alpha, beta, G.X(i,1,1))			*
			pow(1. - G.X(i,1,1) , alpha - 1./2.)				* 
			pow(1. + G.X(i,1,1) , beta - 1./2.)					* 
			G.W(i,1,1);
		}
		integrale /= ( pow((long double)2., alpha + beta + 2));
		return integrale;
	}
}

//integrale di 2 armoniche ipersferiche, all'ordine "ordine", per indice i
inline long double hyper_int (const HyperInt & G,
							  const unsigned short i,
							  const hyper H1, 
							  const hyper H2){
	unsigned short ordine = G.N();
	
	assert (ordine != 0 && i < H1.N() && i < H2.N());
	
	if (H1.l(i) != H2.l(i)) return 0.;
	
	else if (H1.m(i) != H2.m(i)) return 0.;
	
	else if (i == 0) return	int_2j (H1.l(i) , H2.l(i) , 
									H1.m(i) , H2.m(i));
	
	//controllo se i numeri quantici nu al passo precedente è uguale
	
	else if (H1.alpha(i) != H2.alpha(i) || H1.beta(i) != H2.beta(i)) {return 0.;}
	else {			
		long double integrale = int_2j (H1.l(i) , H2.l(i) , 
										H1.m(i) , H2.m(i));
		if ((integrale) != 0) return hyper_int_ripelle(G, i, H1, H2)*integrale;
		else return integrale;
	}
}

//integrale di 2 armoniche ipersferiche, all'ordine "ordine", per tutti gli indici
inline long double hyper_int (const HyperInt & G,
							  const hyper H1, 
							  const hyper H2){
	unsigned short ordine = G.N();
	
	assert (ordine != 0 && H1.N() != 0 && H2.N() != 0 || H1.N() == H2.N());
	
	long double integrale = 1.;
	for (short i = 0; i < H1.N() && integrale != 0.  ; ++i) {
		integrale *= hyper_int ( G, i, H1, H2);
	}
	return integrale;
}


//integrale di 3 funzioni di ripelle, all'ordine "ordine", per indice i
inline long double hyper_int_ripelle (const HyperInt & G,
									  const unsigned short i,
									  const hyper H1, 
									  const hyper H2,
									  const hyper H3){
	unsigned short ordine = G.N();
	
	assert (ordine != 0 && i < H1.N() && i < H2.N() && i < H3.N());
	
	///vedere
	
	if (i == 0) return 1.;
	
	else if (dis_triang(H2.k(i), H1.k(i), H3.k(i)) == false) return 0.;
	
	
	else { 
		long double integrale = 0.;
		unsigned short dub_alpha =  H1.k(i-1) + H2.k(i-1) + H3.k(i-1) + 3*i - 2; //attenzione: sarebbe 3i -5 ma il mio i parte da 0
		unsigned short dub_beta = H1.l(i) + H2.l(i) + H3.l(i) + 1;
		unsigned short n[3] =  { H1.n(i) , H2.n(i) , H3.n(i) };
		long double alpha[3] =  { H1.alpha(i) , H2.alpha(i) , H3.alpha(i) };
		long double beta[3] =  { H1.beta(i) , H2.beta(i) , H3.beta(i) };
		
		if(dub_alpha %2 == 0){ //pari pari -> gauss
			if(dub_beta %2 == 0){
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short I = 0 ; I < ordine; ++I) {
					long double temp = 1.;				
					for (short j = 0 ; j < 3; ++j) {
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}
				integrale /= pow ((long double) 2., (long double) dub_alpha / 2. + dub_beta/2. + 2.);
				return integrale;
			}
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = 1.;
					
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) ,(long double) dub_alpha/2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) ,(long double) dub_beta/2. - 1./2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}
				integrale /= pow ((long double) 2., (long double) dub_alpha / 2. + dub_beta/2. + 2.);
				return integrale;
				
			}
		}
		
		else {
			// alpha semintero, beta intero
			if(dub_beta %2 == 0){  // legendre 2 non modificata
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2 )) *
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2. - 1./2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}
				integrale /= pow ((long double) 2., (long double) dub_alpha / 2. + dub_beta/2. + 2.);
				return integrale;
			}
			// alpha semintero, beta semintero
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *						
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
						
					}
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2. - 1./2.)	*
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2. - 1./2.)	*
								  G.W(i,dub_alpha,dub_beta));
				}
				integrale /= pow ((long double) 2., (long double) dub_alpha / 2. + dub_beta/2. + 2.);
				return integrale;
				
				
			}
			
			
		}
	}
}

//integrale di 3 funzioni di ripelle, all'ordine "ordine", per indice i
inline long double hyper_int (const HyperInt & G,
							  const unsigned short i,
							  const hyper H1, 
							  const hyper H2,
							  const hyper H3){
	unsigned short ordine = G.N();
	assert (ordine != 0 && i < H1.N() && i < H2.N() && i < H3.N());
	
	if (i == 0) return	int_3j (H1.l(i) ,  H2.l(i) , H3.l(i) , 
								H1.m(i) ,  H2.m(i) , H3.m(i));
	
	else {		
		
		long double integrale = int_3j (H1.l(i) ,  H2.l(i) , H3.l(i) , 
										H1.m(i) ,  H2.m(i) , H3.m(i));
		if (integrale != 0) return hyper_int_ripelle(G, i, H1, H2, H3) * integrale;
		else return integrale;
	}
}

//integrale di 3 armoniche ipersferiche, all'ordine "ordine"

inline long double hyper_int_ridotto (const HyperInt & G,
									  const hyper H1, 
									  const hyper H2, 
									  const hyper H3){	
	unsigned short ordine = G.N();
	unsigned short N = H1.N();
	assert (ordine != 0 && N == H2.N() && N == H3.N() && N > 0);
	if ((dis_triang_H_l(H2,H1,H3))){
		
		long double integrale = int_3j_ridotto (H1.l(0) , H2.l(0) , H3.l(0));
		
#pragma omp parallel for reduction(*:integrale) schedule( guided )
		for (short i = 1; i< N ; ++i) {
			long double C = int_3j_ridotto (H1.l(i),  H2.l(i) , H3.l(i));
			if (C != 0) C *= hyper_int_ripelle(G, i, H1, H2, H3);
			else i = N;
			integrale *= C;
		}
		return integrale;
	}
	else return 0.;
}

inline long double hyper_int_ridotto_p (const HyperInt & G,
										const hyper H1, 
										const hyper H2, 
										const hyper H3){	
	unsigned short ordine = G.N();
	unsigned short N = H1.N();
	assert (ordine != 0 && N == H2.N() && N == H3.N() && N > 0);
	if ((dis_triang_H_l(H2,H1,H3))){
		
		long double integrale = 1.;
		
#pragma omp parallel for reduction(*:integrale) schedule( guided )
		for (short i = 1; i< N ; ++i) {
			long double C = hyper_int_ripelle(G, i, H1, H2, H3);
			if (C == 0) i = N;
			integrale *= C;
		}
		return integrale;
	}
	else return 0.;
}

inline long double hyper_int_ridotto_l (const hyper H1, 
										const hyper H2, 
										const hyper H3){	
	unsigned short N = H1.N();
	assert (N == H2.N() && N == H3.N() && N > 0);
	if ((dis_triang_H_l(H2,H1,H3))){
		
		long double integrale = int_3j_ridotto (H1.l(0) , H2.l(0) , H3.l(0));
		
#pragma omp parallel for reduction(*:integrale) schedule( guided )
		for (short i = 1; i< N ; ++i) {
			long double C = int_3j_ridotto (H1.l(i),  H2.l(i) , H3.l(i));
			integrale *= C;
			if (C == 0) i = N;
		}
		return integrale;
	}
	else return 0.;
}


inline long double hyper_int_m (const hyper H1, 
								const hyper H2, 
								const hyper H3){	
	unsigned short N = H1.N();
	assert (N == H2.N() && N == H3.N() && N > 0);
	if (dis_triang_m(H2.m(),H1.m(),H3.m())){
		
		long double integrale = 1.;
		
#pragma omp parallel for reduction(*:integrale) schedule( guided )
		for (short i = 0; i< N ; ++i) {
			integrale *= int_3j_m (H1.l(i) ,  H2.l(i) , H3.l(i) , 
								   H1.m(i) ,  H2.m(i) , H3.m(i));
			if (integrale == 0) i = N;
		}
		return integrale;
	}
	else return 0.;
}

inline long double hyper_int (const HyperInt & G,
							  const hyper H1, 
							  const hyper H2, 
							  const hyper H3){	
	long double integrale = hyper_int_m (H1,H2,H3);
	if (integrale != 0.) integrale *= hyper_int_ridotto (G,H1,H2,H3);
	return integrale;
}

#endif
