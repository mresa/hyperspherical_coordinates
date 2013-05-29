/*
 *  raggio.h
 *  Tesi
 *
 *  Created by Marco Resa on 01/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef RAGGIO_H
#define RAGGIO_H
#include "tipi_definiti.h"
#include "algebra_matrici.h"
#include "interfaces_ublas.h"
#include <boost/numeric/ublas/io.hpp>
#include <vector>

//integrazione di 2 armoniche ipersferiche per l'indice i con una potenza di seni e coseni nel mezzo
inline long double hyper_int_ripelle_trigon (const HyperInt & G,
											 const unsigned short i, // indice della coordinata di jacobi su cui sto integrando. < N
											 const hyper H1, 
											 const hyper H2,
											 const short & coeff_cos,
											 const short & coeff_sin){
	unsigned short ordine = G.N();
	
	assert (ordine != 0 && i < H1.N() && i < H2.N() );
	
	///vedere
	
	if (i == 0) return 1.;//non c'è integrazione che tenga
	
	else { 
		long double integrale = 0.;
		unsigned short dub_alpha =  H1.k(i-1) + H2.k(i-1) + coeff_sin + 3*i - 2; //attenzione: sarebbe 3i -5 ma il mio i parte da 0
		unsigned short dub_beta = H1.l(i) + H2.l(i) + coeff_cos + 1;
		unsigned short n[2] =  { H1.n(i) , H2.n(i)};
		long double alpha[2] =  { H1.alpha(i) , H2.alpha(i) };
		long double beta[2] =  { H1.beta(i) , H2.beta(i) };
		
		if(dub_alpha %2 == 0){ //pari pari -> gauss
			if(dub_beta %2 == 0){
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short I = 0 ; I < ordine; ++I) {
					long double temp = 1.;				
					for (short j = 0 ; j < 2; ++j) {
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}
			}
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = 1.;
					
					for (short j = 0 ; j < 2; ++j) {
						
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) ,(long double) dub_alpha/2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) ,(long double) dub_beta/2. - 1./2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}				
			}
		}
		
		else {
			// alpha semintero, beta intero
			if(dub_beta %2 == 0){  // legendre 2 non modificata
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 2; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2 )) *
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
					}
					
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2. - 1./2.)	* 
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2.)	* 
								  G.W(i,dub_alpha,dub_beta));
				}
			}
			// alpha semintero, beta semintero
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 2; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *						
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
						
					}
					integrale += (temp *
								  pow(1. - G.X(i,dub_alpha,dub_beta) , (long double) dub_alpha/2. - 1./2.)	*
								  pow(1. + G.X(i,dub_alpha,dub_beta) , (long double) dub_beta/2. - 1./2.)	*
								  G.W(i,dub_alpha,dub_beta));
				}
			}			
		}

		integrale /= pow ((long double) 2., (long double) dub_alpha / 2. + dub_beta/2. + 2.);
		return integrale;

	}
}


inline MATRICI_SPARSE rho_quad_m(const unsigned short & N, 
								 const unsigned short & m_max,
								 const Gauss_Laguerre & G, 
								 const long double & beta) {	
	
	MATRICI_SIMMETRICHE_UP V1(ublas::zero_matrix<tipo_matrici>(m_max+1));
	long double alpha = 3.*N - 1.;
	
	//integrazioni sulla base di stati
	for (int k = 0; k < G.N() ; ++k){
		std::vector <long double> L;
		for (int i=0; i< m_max + 1 ; ++i)	{
			L.push_back( laguerre ( i , alpha , G.X(k) ) );
		}
		
		long double ro_quad = pow( G.X(k)/beta, 2);	
		
		long double C = G.W(k) *  pow( G.X(k) , (alpha - G.ALPHA()) );

		
		for (int i=0; i< m_max + 1 ; ++i)	{
			for (int j=i; j< m_max + 1; ++j)	{
				long double C1 = sqrt( boost::math::tgamma_delta_ratio(i + 1, alpha) * boost::math::tgamma_delta_ratio(j + 1, alpha) );
				V1(i,j) += C * C1 * L[j] * L[i] * ro_quad ;
		}
		}
	}
	MATRICI_SPARSE RISULTATO(m_max+1,m_max+1);
	for (int i=0; i< m_max + 1 ; ++i)	{
		for (int j=i; j< m_max + 1; ++j)	{
			if (V1(i,j)) {
				RISULTATO(i,j) = V1(i,j);
				RISULTATO(j,i) = V1(i,j);
			}	
		}
	}
	

// controllo se definita positiva
//	MATRICI V2(V1);
//	for (unsigned short i = 1; i < V1.size1() +1 ; ++i) {
//		ublas::matrix_range<ublas::matrix<double> > mr (V2, ublas::range (0, i), ublas::range (0, i));
//		ublas::matrix<double> V3(mr);
//		std::cerr << "determinanti dei minori" <<  ublas::lu_det (V3) <<std::endl;
//	}
//	std::cerr << "rho " <<  V1 <<std::endl;

	return RISULTATO;
}

inline MATRICI_SPARSE x_i_quad_k(const Accoppiamenti_Iperangolari_K_max & K,
								 const unsigned short & indice, // della coordinata dijacobi
								 const HyperInt & G,
								 const unsigned short & dimensione_calcolata) {
	unsigned short N = K.N();
	
	assert(indice < N);
	
	MATRICI_SPARSE M(K.D(),K.D());
	
	for(short i = 0 ; i < K.D() ; ++i){
		hyper H1_0 (K.n(i),K.l(i),K.m(i,0));
		for(short j = i ; j < K.D() ; ++j){
			if( i >= dimensione_calcolata || j >= dimensione_calcolata){
				if(delta_vec(K.l(i),K.l(j)) && delta_vec_fino(K.n(i),K.n(j),indice)){ //è diagonale sulle l
					hyper H2_0 (K.n(j),K.l(j),K.m(j,0));
					long double integrale = hyper_int_ripelle_trigon(G, indice , H2_0, H1_0, 2,0); //cos^2
					if (integrale != 0) {
						for (short ii = indice + 1; ii< N ; ++ii) {
							integrale *= hyper_int_ripelle_trigon(G, ii,H2_0,H1_0,0,2); //sin^2
							if (integrale == 0) ii = N;
						}
					}
					if(integrale != 0){
						long double temp = 0.;
						for(short d1 = 0; d1 < K.DM(i); ++d1){
							long double temp1 = 0.;
							for(short d2 = 0; d2 < K.DM(j); ++d2){
								temp1 += K.cm(j,d2)*delta_vec(K.m(i,d1),K.m(j,d2));
							}
							temp += K.cm(i,d1) * temp1;
						}
						integrale*=temp;
						if(integrale){
							M(i,j) = integrale;
							M(j,i) = integrale;
						}
					}
				}
			}
		}
	}
	//	std::cerr << "indice = " << indice << "M _ x2" <<  M <<std::endl;
	//controllo se definita positiva
	//	MATRICI V2(M);
	//	for (unsigned short i = 1; i < V2.size1() +1 ; ++i) {
	//		ublas::matrix_range<ublas::matrix<double> > mr (V2, ublas::range (0, i), ublas::range (0, i));
	//		ublas::matrix<double> V3(mr);
	//		std::cerr << "determinanti dei minori" <<  ublas::lu_det (V3) <<std::endl;
	//	}
	return M;
}

inline MATRICI_SPARSE x_i_quad_k(const Accoppiamenti_Iperangolari_K_max & K,
								 const unsigned short & indice, // della coordinata dijacobi
								 const HyperInt & G) {
	return x_i_quad_k(K,indice,G,0);
}


inline MATRICI_SPARSE x_i_j_k(const Accoppiamenti_Iperangolari_K_max & K,
							  const unsigned short & indice, // della coordinata dijacobi
							  const unsigned short & jndice, // della coordinata dijacobi MAGGIORE!!
							  const HyperInt & G,
							  const unsigned short & dimensione_calcolata) {
	unsigned short N = K.N();
	unsigned short indice_t; //minore
	unsigned short jndice_t; //maggiore
	assert(indice < N);
	assert(jndice < N);
	
	if (indice < jndice){
		indice_t = indice;
		jndice_t = jndice;
	}
	
	else if(indice > jndice){
		indice_t = jndice;
		jndice_t = indice;
	}
	
	else {
		std::cerr << "chiamato x_ij con indici uguali";
		return x_i_quad_k(K,indice,G,dimensione_calcolata);
	}
	
	MATRICI_SPARSE M(K.D(),K.D());
	
	for(unsigned short i = 0 ; i < K.D() ; ++i){
		hyper H1_0 (K.n(i),K.l(i),K.m(i,0));
		for(short j = i ; j < K.D() ; ++j){
			if( i >= dimensione_calcolata || j >= dimensione_calcolata){
				if(delta_vec_tranne(K.l(i),K.l(j),indice_t,jndice_t) && delta_vec_fino(K.n(i),K.n(j),indice)){ //è diagonale sulle l,tranne in i e j e devono essere uguali gli indici fino al desiderato
					hyper H2_0 (K.n(j),K.l(j),K.m(j,0));
					long double integrale = hyper_int_ripelle_trigon(G, indice_t,H2_0,H1_0,1,0);
					if (integrale != 0) {
						for (short ii = indice_t + 1; ii< jndice_t ; ++ii) {
							integrale *= hyper_int_ripelle_trigon(G, ii,H2_0,H1_0,0,1);
							if (integrale == 0) ii = N;
						}
					}
					if (integrale != 0) {
						integrale *= hyper_int_ripelle_trigon(G, jndice_t,H2_0,H1_0,1,1);
						for (short ii = jndice_t + 1; ii< N ; ++ii) {
							integrale *= hyper_int_ripelle_trigon(G, ii,H2_0,H1_0,0,2);
							if (integrale == 0) ii = N;
						}
					}
					if(integrale != 0){
						long double temp = 0.;
						for(short d1 = 0; d1 < K.DM(i); ++d1){
							long double temp1 = 0.;
							for(short d2 = 0; d2 < K.DM(j); ++d2){
								if(delta_vec_tranne(K.m(i,d1),K.m(j,d2),indice_t,jndice_t)){
									long double temp2 = int_3j(K.l(j)[indice_t],1,K.l(i)[indice_t],K.m(j,d2)[indice_t],0,K.m(i,d1)[indice_t])*
									int_3j(K.l(j)[jndice_t],1,K.l(i)[jndice_t],K.m(j,d2)[jndice_t],0,K.m(i,d1)[jndice_t]);
									
									temp2 -= int_3j(K.l(j)[indice_t],1,K.l(i)[indice_t],K.m(j,d2)[indice_t],1,K.m(i,d1)[indice_t])*
									int_3j(K.l(j)[jndice_t],1,K.l(i)[jndice_t],K.m(j,d2)[jndice_t],-1,K.m(i,d1)[jndice_t]);
									
									temp2 -= int_3j(K.l(j)[indice_t],1,K.l(i)[indice_t],K.m(j,d2)[indice_t],-1,K.m(i,d1)[indice_t])*
									int_3j(K.l(j)[jndice_t],1,K.l(i)[jndice_t],K.m(j,d2)[jndice_t],1,K.m(i,d1)[jndice_t]);

									temp2 *= 4*pi/3;
									temp1 += K.cm(j,d2)*temp2;
								}
							}
							temp += K.cm(i,d1) * temp1;
						}
						integrale*=temp;
						if (integrale) {
							M(i,j) = integrale;
							M(j,i) = integrale;
						}
					}
				}
			}
		}
	}
	//	std::cerr << "indice = " << indice << "   jndice = " << jndice << " M _ x2" <<  M <<std::endl;
	//controllo se definita positiva
	//	MATRICI V2(M);
	//	for (unsigned short i = 1; i < V2.size1() +1 ; ++i) {
	//		ublas::matrix_range<ublas::matrix<double> > mr (V2, ublas::range (0, i), ublas::range (0, i));
	//		ublas::matrix<double> V3(mr);
	//		std::cerr << "determinanti dei minori" <<  ublas::lu_det (V3) <<std::endl;
	//	}
	return M;
}

inline MATRICI_SPARSE x_i_j_k(const Accoppiamenti_Iperangolari_K_max & K,
							  const unsigned short & indice, // della coordinata dijacobi
							  const unsigned short & jndice, // della coordinata dijacobi MAGGIORE!!
							  const HyperInt & G) {
	return x_i_j_k(K,indice,jndice,G,0);
}

inline V_MATRICI_SPARSE x_quad_k(const Accoppiamenti_Iperangolari_K_max & K,
								 const HyperInt & G,
								 const unsigned short & dimensione_calcolata){
	V_MATRICI_SPARSE x_quad_k_t;
	for (short i = 0; i<K.N(); ++i) {
		x_quad_k_t.push_back(x_i_quad_k(K,i,G,dimensione_calcolata));
	}
	return x_quad_k_t;
}

inline V_MATRICI_SPARSE x_quad_k(const Accoppiamenti_Iperangolari_K_max & K,
								 const HyperInt & G){
	return x_quad_k(K,G,0);
}


inline V_MATRICI_SPARSE x_ij_k(const Accoppiamenti_Iperangolari_K_max & K,
							   const HyperInt & G,
							   const unsigned short & dimensione_calcolata){
	V_MATRICI_SPARSE x_ij_k_t;
	for (short i = 0; i<K.N(); ++i) {
		for (short j = i+1; j<K.N(); ++j) {
			x_ij_k_t.push_back(x_i_j_k(K,i,j,G,dimensione_calcolata));
		}
	}
	return x_ij_k_t;
}


inline V_MATRICI_SPARSE x_ij_k(const Accoppiamenti_Iperangolari_K_max & K,
							   const HyperInt & G){
	return x_ij_k(K,G,0);
}

//costruisco il vettore r_i_quad e la matrice r_ij_quad

inline VETTORI r_quadro(const V_V_V_MATRICI_SPARSE & R_QUADRO_E,
						const VETTORI & v){
	VETTORI r_quad_t(R_QUADRO_E.size());
	std::cerr << "R_QUADRO_E.size() = " << R_QUADRO_E.size() << std::endl;
	for(short i = 0 ; i < R_QUADRO_E.size() ; ++i){
		r_quad_t(i) = ublas::inner_prod(v, v_mult(R_QUADRO_E[i],v));
	}
	return ublas::sqrt_vec(r_quad_t);
}

inline MATRICI_SIMMETRICHE_UP r_delta_quadro(const V_V_V_MATRICI_SPARSE & R_D_QUADRO_E,
											 const unsigned short & A_t,
											 const VETTORI & v){
	MATRICI_SIMMETRICHE_UP r_delta_quad_t(A_t,A_t);
	unsigned short n = 0;
	for (short i = 0; i < A_t; ++i) {
		for (short j = i; j < A_t; ++j) {
			VETTORI v1(ublas::v_mult(R_D_QUADRO_E[n],v));
			r_delta_quad_t(i,j) = ublas::inner_prod(v, v1);
			++n;
		}
	}
	return ublas::sqrt_mat_sim_up(r_delta_quad_t);
}

inline tipo_matrici valor_medio(const V_V_MATRICI_SPARSE & M,
								const VETTORI & v){

	if (M.size()) {
		return ublas::inner_prod(v, v_mult(M,v));
	}
	else {
		return 0.;
	}

}


#endif
