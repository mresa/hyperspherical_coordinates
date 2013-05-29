/*
 *  potenziale
 *  Tesi
 *
 *  Created by Marco Resa on 29/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef POTENZIALE_H
#define POTENZIALE_H

#include "interazione.h"
#include "tipi_definiti.h"
#include "hyperint.h"
#include <vector>
#include <omp.h> //calcolo parallelo OMP


void potenziale(const unsigned short & N, 
				const unsigned short & m_max, 
				const Gauss_Laguerre & G , 
				const std::vector <long double> & V,
				MATRICI_SPARSE * V1) {
	assert( (*V1).size1() == m_max + 1);
	assert( (*V1).size2() == m_max + 1);
	
	MATRICI_SIMMETRICHE_UP V1_t(ublas::zero_matrix<tipo_matrici>((*V1).size1()));
	long double alpha = 3.*N - 1.;
	
	for (int k = 0; k < G.N() ; ++k)	{
		
		std::vector <long double> L;
		
		//costruisco i laguerre che userò
		
		for (int i=0; i< m_max + 1 ; ++i)	{
			L.push_back( laguerre ( i , alpha , G.X(k) ) );
		}
		
		long double C = G.W(k) *  pow( G.X(k) , (alpha - G.ALPHA()) );
		for (int i=0; i< m_max + 1 ; ++i)	{
			for (int j=i; j< m_max + 1; ++j)	{
				long double C1 = sqrt( boost::math::tgamma_delta_ratio(i + 1, alpha) * boost::math::tgamma_delta_ratio(j + 1, alpha) );
				V1_t(i,j) += C * C1 * L[j] * L[i] * V[k] ;
			}
		}
	}
	for (int i=0; i< m_max + 1 ; ++i)	{
		for (int j=i; j< m_max + 1; ++j)	{
			if (V1_t(i,j)) {
				(*V1)(i,j) = V1_t(i,j);
				(*V1)(j,i) = V1_t(i,j);
			}	
		}
	}
	
}


inline std::vector < long double > V_n_Gauss_Laguerre (const unsigned short & N,
													   const HyperInt & HI,
													   const Gauss_Laguerre & GL,
													   const unsigned short & n,
													   const long double & BETA,
													   const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione(ro));
		}
	}
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){

			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
					for (short i = 0 ; i < ordine; ++i) {
						long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
						jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));

						integrale += (temp *
									  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
									  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
									  HI.W(i,dubalpha,dubbeta));
					}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	return temp;
}

void potenziale (const unsigned short & N,
				 const HyperInt & HI,
				 const Gauss_Laguerre & G,
				 const unsigned short & m_max,
				 const unsigned short & n,
				 const long double & beta,
				 const long double & C_IJ,
				 MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}


inline std::vector < long double > V_n_C_Gauss_Laguerre (const unsigned short & N,
														 const HyperInt & HI,
														 const Gauss_Laguerre & GL,
														 const unsigned short & n,
														 const long double & BETA,
														 const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_coulombiana(ro));
		}
	}
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_coulombiana(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_coulombiana(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	/*else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_coulombiana(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	return temp;
}


void potenziale_coulombiano (const unsigned short & N,
							 const HyperInt & HI,
							 const Gauss_Laguerre & G,
							 const unsigned short & m_max,
							 const unsigned short & n,
							 const long double & beta,
							 const long double & C_IJ,
							 MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_C_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}


V_MATRICI_SPARSE potenziale (const unsigned short & N,
							 const HyperInt & HI,
							 const Gauss_Laguerre & G,
							 const unsigned short & m_max,
							 const unsigned short & n,
							 const long double & beta,
							 const std::vector < long double > & C_IJ) {
	V_MATRICI_SPARSE V;
	for(short i = 0 ; i < C_IJ.size() ; ++i) {
		MATRICI_SPARSE V_temp (m_max + 1,m_max + 1);					
		potenziale(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp);
		V.push_back(V_temp);
	}
	return V;
}

V_MATRICI_SPARSE potenziale_coulombiano (const unsigned short & N,
										 const HyperInt & HI,
										 const Gauss_Laguerre & G,
										 const unsigned short & m_max,
										 const unsigned short & n,
										 const long double & beta,
										 const std::vector < long double > & C_IJ,
										 const std::vector < bool > & cariche_ij) {
	V_MATRICI_SPARSE V;
	for(short i = 0 ; i < C_IJ.size() ; ++i) {
		if (cariche_ij[i]) {
			MATRICI_SPARSE V_temp (m_max + 1,m_max + 1);					
			potenziale_coulombiano(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp);
			V.push_back(V_temp);
		}
		else { //non lo calcolo
			MATRICI_SPARSE V_temp (m_max + 1,m_max + 1);
			V.push_back(V_temp);
		}
	}
	return V;
}

V_MATRICI_SPARSE potenziale_coulombiano (const unsigned short & N,
										 const HyperInt & HI,
										 const Gauss_Laguerre & G,
										 const unsigned short & m_max,
										 const unsigned short & n,
										 const long double & beta,
										 const std::vector < long double > & C_IJ) {
	V_MATRICI_SPARSE V;
	for(short i = 0 ; i < C_IJ.size() ; ++i) {
		MATRICI_SPARSE V_temp (m_max + 1,m_max + 1);					
		potenziale_coulombiano(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp);
		V.push_back(V_temp);
	}
	return V;
}

//spinoriale

inline std::vector < long double > V_n_s_0_Gauss_Laguerre (const unsigned short & N,
														   const HyperInt & HI,
														   const Gauss_Laguerre & GL,
														   const unsigned short & n,
														   const long double & BETA,
														   const long double & C_IJ) {
	std::vector <long double> temp;
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_s_0(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_s_0(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_s_0(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_s_0(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	return temp;
}

inline std::vector < long double > V_n_s_1_Gauss_Laguerre (const unsigned short & N,
														   const HyperInt & HI,
														   const Gauss_Laguerre & GL,
														   const unsigned short & n,
														   const long double & BETA,
														   const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_s_1(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_s_1(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_s_1(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_s_1(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	return temp;
}

void potenziale_s_0 (const unsigned short & N,
					 const HyperInt & HI,
					 const Gauss_Laguerre & G,
					 const unsigned short & m_max,
					 const unsigned short & n,
					 const long double & beta,
					 const long double & C_IJ,
					 MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_s_0_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}

void potenziale_s_1 (const unsigned short & N,
					 const HyperInt & HI,
					 const Gauss_Laguerre & G,
					 const unsigned short & m_max,
					 const unsigned short & n,
					 const long double & beta,
					 const long double & C_IJ,
					 MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_s_1_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}



V_V_MATRICI_SPARSE potenziale_s (const unsigned short & N,
								 const HyperInt & HI,
								 const Gauss_Laguerre & G,
								 const unsigned short & m_max,
								 const unsigned short & n,
								 const long double & beta,
								 const std::vector < long double > & C_IJ) {
	V_V_MATRICI_SPARSE V;
	V_MATRICI_SPARSE V_s_0;
	V_MATRICI_SPARSE V_s_1;
	for(short i = 0 ; i < C_IJ.size() ; ++i) {
		MATRICI_SPARSE V_temp_s_0 (m_max + 1,m_max + 1);					
		MATRICI_SPARSE V_temp_s_1 (m_max + 1,m_max + 1);					
		potenziale_s_0(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_s_0);
		potenziale_s_1(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_s_1);
		V_s_0.push_back(V_temp_s_0);
		V_s_1.push_back(V_temp_s_1);
	}
	V.push_back(V_s_0);
	V.push_back(V_s_1);
	return V;
}


//isospinoriale

inline std::vector < long double > V_n_st_00_Gauss_Laguerre (const unsigned short & N,
															 const HyperInt & HI,
															 const Gauss_Laguerre & GL,
															 const unsigned short & n,
															 const long double & BETA,
															 const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_st_00(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_st_00(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_00(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_00(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 	
	return temp;
}

inline std::vector < long double > V_n_st_01_Gauss_Laguerre (const unsigned short & N,
															 const HyperInt & HI,
															 const Gauss_Laguerre & GL,
															 const unsigned short & n,
															 const long double & BETA,
															 const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_st_01(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_st_01(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_01(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_01(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	
	return temp;
}

inline std::vector < long double > V_n_st_10_Gauss_Laguerre (const unsigned short & N,
															 const HyperInt & HI,
															 const Gauss_Laguerre & GL,
															 const unsigned short & n,
															 const long double & BETA,
															 const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_st_10(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_st_10(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_10(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_10(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 
	
	return temp;
}

inline std::vector < long double > V_n_st_11_Gauss_Laguerre (const unsigned short & N,
															 const HyperInt & HI,
															 const Gauss_Laguerre & GL,
															 const unsigned short & n,
															 const long double & BETA,
															 const long double & C_IJ) {
	std::vector <long double> temp;
	
	if (N==1){
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			temp.push_back(interazione_st_11(ro));
		}
	}
/*	else {
		unsigned short l = 0;
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = l+1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.;
		
		long double jac_alpha = alpha;
		long double jac_beta = beta + l/2.;
		
		Integrale_Pippo G(HI,dubalpha,dubbeta);
		long double C = potential_basis_normalizzazione(N,n,l)/(hyper0(N-1)*pow((long double)2., alpha + beta + 2 ));
		C *= 2. * sqrt(pi); //dovuto alla parte centrale : è l'armonica sferica di grado 0
		for (short k = 0 ; k < GL.N() ; ++k){
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
#pragma omp parallel for reduction(+:integrale) schedule( guided )
			for (short i = 0 ; i < G.N(); ++i) {
				integrale += (G.W(i)											*
							  jacobi( n, jac_alpha, jac_beta, G.X(i) )			*
							  interazione_st_11(ro*sqrt((1. + G.X(i))/2.)));
			}
			integrale *= C;
			temp.push_back(integrale);
		}
	}*/
	else {
		unsigned short ordine = HI.N();
		unsigned short dubalpha = (3*N - 5);
		unsigned short dubbeta = 1;
		long double alpha = dubalpha/2.;
		long double beta = dubbeta/2.; //sempre dispari!
		
		long double jac_alpha = alpha;
		long double jac_beta = beta;
		
		long double C = 2. * sqrt(pi)/(hyper0(N-1)*sqrt(pow((long double)2., 3*N )));
		for (short k = 0 ; k < GL.N() ; ++k){
			
			long double ro = GL.X(k)*C_IJ/BETA;
			long double integrale = 0.;
			if(dubalpha %2 == 0){ 
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_11(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) ,dubalpha))	* 
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) ,dubbeta - 1))	* 
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
				
			}
			
			else {
				// alpha semintero, beta semintero
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = sqrt(pow ((long double) 2., jac_alpha + jac_beta + 2.)) *						
					jacobi_N( n, jac_alpha, jac_beta, HI.X(i,dubalpha,dubbeta))*interazione_st_11(ro*sqrt((1. + HI.X(i,dubalpha,dubbeta))/2.));
					
					integrale += (temp *
								  sqrt(pow(1. - HI.X(i,dubalpha,dubbeta) , dubalpha - 1))	*
								  sqrt(pow(1. + HI.X(i,dubalpha,dubbeta) , dubbeta - 1))	*
								  HI.W(i,dubalpha,dubbeta));
				}
				integrale *= C;
				temp.push_back(integrale);
			}				
		}			
	} 	
	return temp;
}

void potenziale_st_00 (const unsigned short & N,
					   const HyperInt & HI,
					   const Gauss_Laguerre & G,
					   const unsigned short & m_max,
					   const unsigned short & n,
					   const long double & beta,
					   const long double & C_IJ,
					   MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_st_00_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}

void potenziale_st_01 (const unsigned short & N,
					   const HyperInt & HI,
					   const Gauss_Laguerre & G,
					   const unsigned short & m_max,
					   const unsigned short & n,
					   const long double & beta,
					   const long double & C_IJ,
					   MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_st_01_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}

void potenziale_st_10 (const unsigned short & N,
					   const HyperInt & HI,
					   const Gauss_Laguerre & G,
					   const unsigned short & m_max,
					   const unsigned short & n,
					   const long double & beta,
					   const long double & C_IJ,
					   MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_st_10_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}

void potenziale_st_11 (const unsigned short & N,
					   const HyperInt & HI,
					   const Gauss_Laguerre & G,
					   const unsigned short & m_max,
					   const unsigned short & n,
					   const long double & beta,
					   const long double & C_IJ,
					   MATRICI_SPARSE * V1) {
	potenziale( N, m_max,  G , V_n_st_11_Gauss_Laguerre(N,HI,G,n,beta,C_IJ), V1 );
}


V_V_MATRICI_SPARSE potenziale_st (const unsigned short & N,
								  const HyperInt & HI,
								  const Gauss_Laguerre & G,
								  const unsigned short & m_max,
								  const unsigned short & n,
								  const long double & beta,
								  const std::vector < long double > & C_IJ) {
	V_V_MATRICI_SPARSE V;
	V_MATRICI_SPARSE V_st_00;
	V_MATRICI_SPARSE V_st_01;
	V_MATRICI_SPARSE V_st_10;
	V_MATRICI_SPARSE V_st_11;
	for(short i = 0 ; i < C_IJ.size() ; ++i) {
		MATRICI_SPARSE V_temp_st_00 (m_max + 1,m_max + 1);					
		MATRICI_SPARSE V_temp_st_10 (m_max + 1,m_max + 1);					
		MATRICI_SPARSE V_temp_st_01 (m_max + 1,m_max + 1);					
		MATRICI_SPARSE V_temp_st_11 (m_max + 1,m_max + 1);					
		potenziale_st_00(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_st_00);
		potenziale_st_01(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_st_01);
		potenziale_st_10(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_st_10);
		potenziale_st_11(N,HI,G,m_max,n,beta,C_IJ[i],&V_temp_st_11);
		
		V_st_00.push_back(V_temp_st_00);
		V_st_01.push_back(V_temp_st_01);
		V_st_10.push_back(V_temp_st_10);
		V_st_11.push_back(V_temp_st_11);
	}
	V.push_back(V_st_00);
	V.push_back(V_st_10);
	V.push_back(V_st_01);
	V.push_back(V_st_11);
	return V;
}


#endif

