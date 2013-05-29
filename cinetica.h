/*
 *  cinetica_ublas.h
 *  Tesi
 *
 *  Created by Marco Resa on 29/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CINETICA_H
#define CINETICA_H


#include <assert.h>
#include <vector>
#include "gauss-laguerre.h"
#include "algebra.h"
#include "tipi_definiti.h"
#include "algebra_matrici.h"
#include "rotazione_iperangolare.h"
#include "costanti_fisiche.h"
#include <omp.h> //calcolo parallelo OMP

void cinetica(const unsigned short & N,
			  const unsigned short & m_max,
			  const Gauss_Laguerre & G,
			  const long double& beta, 
			  const long double& MASSA_RIFERIMENTO,
			  MATRICI_SPARSE * T1,
			  MATRICI_SPARSE * T2) {
	
	assert(N>0);
	assert((*T1).size1() == m_max + 1);
	assert((*T1).size2() == m_max + 1);
	assert((*T2).size1() == m_max + 1);
	assert((*T2).size2() == m_max + 1);
	long double B = 41.47;//pow(h_bar_c, 2)/MASSA_RIFERIMENTO;
	B *= - pow(beta, 2);
	
	unsigned short alpha = 3*N - 1;
	//memorizzo tutti i laguerre che userò in una matrice L[i][k]
	std::vector <vettore_ldouble> L;
	std::vector <long double> P1;
	std::vector <long double> P2;
	std::vector <long double> P3;
	for (int k = 0; k < G.N() ; ++k)	{
		std::vector <long double> L1;
		for (int i=0; i< (*T1).size1() ; ++i)	{
			L1.push_back( laguerre ( i , alpha , G.X(k) ) );
		}
		L.push_back(L1);
		P1.push_back(pow( G.X(k) , 2 ));
		P2.push_back(pow( G.X(k) , (alpha - G.ALPHA()) ));
		P3.push_back(pow( G.X(k) , (alpha - G.ALPHA() - 2 ) ));		
	}
	
	for (int i=0; i < m_max + 1 ; ++i)	{		
		
		long double C1 = sqrt( boost::math::tgamma_delta_ratio(i + 1, alpha));
		
		for (int j=i; j< m_max + 1; ++j)	{
			long double C2 = sqrt( boost::math::tgamma_delta_ratio(j + 1, alpha) );
			
			long double integrale_1 = 0.;
			
#pragma omp parallel for reduction(+:integrale_1) schedule( guided )
			for (int k = 0; k < G.N() ; ++k){
				long double condizionale;
				if (i==0) condizionale = 0.;
				else condizionale = ( i + alpha) / P1[k] * L[k][i-1];
				
				integrale_1 += G.W(k) * L[k][j] *  P2[k] * ( ( - ( alpha + 2.*i)/( 2. * G.X(k) ) - i / P1[k] )* L[k][i] + condizionale ); 
			}
			
			long double integrale_2 = 0.;
			
#pragma omp parallel for reduction(+:integrale_2) schedule( guided )
			for (int k = 0; k < G.N() ; ++k){
				
				integrale_2 += G.W(k) * P3[k] * L[k][j] * L[k][i]; 
			}
			
			if(integrale_1 || i==j){
				integrale_1 *= C1 * C2;
				integrale_1 += Kdelta(i, j) / 4.;
				(*T1)(i,j) = integrale_1;
				(*T1)(j,i) = integrale_1;
			}
		
			if(integrale_2){
				integrale_2 *= C1 * C2;
				(*T2)(i,j) = integrale_2;
				(*T2)(j,i) = integrale_2;
			}
		}
	}
	(*T1) *= B;
	(*T2) *= B;
}

void nabla1(const unsigned short & N,
			const unsigned short & m_max,
			const Gauss_Laguerre & G,
			const long double& beta, 
			const long double& MASSA_RIFERIMENTO,
			MATRICI_SPARSE * T3) {
	
	assert(N>0);
	assert((*T3).size1() == m_max + 1);
	assert((*T3).size2() == m_max + 1);
	long double B = pow(h_bar_c, 2)/MASSA_RIFERIMENTO;
	B *= - pow(beta, 2);
	unsigned short alpha = 3*N - 1;
	//memorizzo tutti i laguerre che userò in una matrice L[i][k]
	std::vector <vettore_ldouble> L;
	std::vector <long double> P1;
	std::vector <long double> P2;
	std::vector <long double> P3;
	for (int k = 0; k < G.N() ; ++k)	{
		std::vector <long double> L1;
		for (int i=0; i< (*T3).size1() ; ++i)	{
			L1.push_back( laguerre ( i , alpha , G.X(k) ) );
		}
		L.push_back(L1);
		P1.push_back(pow( G.X(k) , 2 ));
		P2.push_back(pow( G.X(k) , (alpha - G.ALPHA()) ));
		P3.push_back(pow( G.X(k) , (alpha - G.ALPHA() - 2 ) ));		
	}
	
	for (int i=0; i< m_max + 1 ; ++i)	{		
		
		long double C1 = sqrt( boost::math::tgamma_delta_ratio(i + 1, alpha));
		
		for (int j=i; j< m_max + 1; ++j)	{
			long double C2 = sqrt( boost::math::tgamma_delta_ratio(j + 1, alpha) );
			
			long double integrale_1 = 0.;
			
#pragma omp parallel for reduction(+:integrale_1) schedule( guided )
			for (int k = 0; k < G.N() ; ++k){
				long double condizionale;
				if (i==0) condizionale = 0.;
				else condizionale = ( i + alpha) / P1[k] * L[k][i-1];
				
				integrale_1 += G.W(k) * L[k][j] *  P2[k] * ( (  ( i + 1)/P1[k]  - 1/( 2. * G.X(k) ) ) * L[k][i] - condizionale ); 
			}
						
			if(integrale_1){
				integrale_1 *= C1 * C2;
				(*T3)(i,j) = integrale_1;
				(*T3)(j,i) = integrale_1;
			}
		}
	}
	(*T3) *= B;
}

V_MATRICI_SPARSE cinetica(const unsigned short & N,
						  const unsigned short & m_max,
						  const Gauss_Laguerre & G,
						  const long double& beta, 
						  const long double& MASSA_RIFERIMENTO) {
	
	assert(N>0);
	V_MATRICI_SPARSE TT;
	MATRICI_SPARSE T1(m_max + 1,m_max + 1);
	MATRICI_SPARSE T2(m_max + 1,m_max + 1);
	
	long double B = pow(h_bar_c, 2)/MASSA_RIFERIMENTO;
	B *= - pow(beta, 2);
	
	unsigned short alpha = 3*N - 1;
	//memorizzo tutti i laguerre che userò in una matrice L[i][k]
	std::vector <vettore_ldouble> L;
	std::vector <long double> P1;
	std::vector <long double> P2;
	std::vector <long double> P3;
	for (int k = 0; k < G.N() ; ++k)	{
		std::vector <long double> L1;
		for (int i=0; i< T1.size1() ; ++i)	{
			L1.push_back( laguerre ( i , alpha , G.X(k) ) );
		}
		L.push_back(L1);
		P1.push_back(pow( G.X(k) , 2 ));
		P2.push_back(pow( G.X(k) , (alpha - G.ALPHA()) ));
		P3.push_back(pow( G.X(k) , (alpha - G.ALPHA() - 2 ) ));		
	}
	
	for (int i=0; i< m_max + 1 ; ++i)	{		
		
		long double C1 = sqrt( boost::math::tgamma_delta_ratio(i + 1, alpha));
		
		for (int j=i; j< m_max + 1; ++j)	{
			long double C2 = sqrt( boost::math::tgamma_delta_ratio(j + 1, alpha) );
			
			long double integrale_1 = 0.;
			
#pragma omp parallel for reduction(+:integrale_1) schedule( guided )
			for (int k = 0; k < G.N() ; ++k){
				long double condizionale;
				if (i==0) condizionale = 0.;
				else condizionale = ( i + alpha) / P1[k] * L[k][i-1];
				
				integrale_1 += G.W(k) * L[k][j] *  P2[k] * ( ( - ( alpha + 2.*i)/( 2. * G.X(k) ) - i / P1[k] )* L[k][i] + condizionale ); 
			}
			
			long double integrale_2 = 0.;
			
#pragma omp parallel for reduction(+:integrale_2) schedule( guided )
			for (int k = 0; k < G.N() ; ++k){
				
				integrale_2 += G.W(k) * P3[k] * L[k][j] * L[k][i]; 
			}
			
			if(integrale_1 || i==j){
				integrale_1 *= C1 * C2;
				integrale_1 += Kdelta(i, j) / 4.;
				T1(i,j) = integrale_1;
				T1(j,i) = integrale_1;
			}
			
			if(integrale_2){
				integrale_2 *= C1 * C2;
				T2(i,j) = integrale_2;
				T2(j,i) = integrale_2;
			}
		}
	}
	T1 *= B;
	T2 *= B;
	TT.push_back(T1);
	TT.push_back(T2);
	return TT;
}

#endif

