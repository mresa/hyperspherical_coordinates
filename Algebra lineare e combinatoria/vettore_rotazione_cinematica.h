
/*
 *  vettore_rotazione_cinematica.h
 *  Tesi
 *
 INPUT : MATRICE DI JACOBI, coppia i j
 
 OUTPUT vettore φ_2, …, φ_N
 
 NOTA BENE in input utilizza l'ordinamento diretto (inserire matrice diretta)
 
 RESTITUISCE ORDINAMENTO INVERSO
 
 *  Created by Marco Resa on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ROT_CIN_H
#define ROT_CIN_H

#include <assert.h>
#include <vector>
#include "tipi_definiti.h"
#include "err.h"
#include "algebra.h"
#include "algebra_matrici.h"
#include "matrice_trasformazione_jacobi.h"




//trovo r_ij = c[1] x_1 + + c[n] x_n
inline std::vector< long double > coefficienti_diretti_r_ij(const ublas::matrix<long double> & M_inversa, 
															const unsigned short& i, 
															const unsigned short& j){
	unsigned short a = M_inversa.size1();
	std::vector< long double > coefficienti_dir;
	for(short l = 0 ; l < a ; ++l) {
		long double temp = (M_inversa(j,l) - M_inversa(i,l));
		coefficienti_dir.push_back(temp);
	}
	assert (std::abs(coefficienti_dir[ a - 1 ])  < err_coordinate );//== 0.); //controllo che il A-esimo sia nullo
	coefficienti_dir.pop_back();	//e lo cancello: ora ne ho n
	return 	coefficienti_dir;
}


inline long double C_rot_cin(const std::vector< long double > & coefficienti_diretti) {
	long double C = 0;
#pragma omp parallel for reduction(+:C) schedule( guided )
	for (short i = 0 ; i < coefficienti_diretti.size()  ; ++i){
		C += pow(coefficienti_diretti[i], 2);
	}
	C = sqrt(C);
	return C;
}

//vettore di dimensione (A-1)*A / 2
//nb per ogni coppia ij è pari a sqrt(m(m_i + m_j)/(2*m_i*m_j)
inline std::vector<long double> C_rot_cin(const ublas::matrix<long double> & M_inversa) {
	unsigned short A = M_inversa.size1();
	
	std::vector<long double> C;
	
	for (short i = 0; i < A ; ++i){
		for (short j = i+1; j < A; ++j){
			C.push_back(C_rot_cin(coefficienti_diretti_r_ij(M_inversa,i,j))) ;
		}
	}
	unsigned short S = ((A-1)*A)/2;
	assert (C.size() == S);
	return C;
}


/*
 *	calcola gli angoli fi_ij del vettore di rotazione cinematica, data la matrice di trasformazione, i, j (vanno da 0 a A-1)
 *
 */


inline std::vector< long double > rot_cin(const ublas::matrix<long double> & M_inversa, 
										  const unsigned short& i, 
										  const unsigned short& j) {
	
	std::vector< long double > coefficienti_dir(coefficienti_diretti_r_ij(M_inversa, i, j));
	unsigned short n = coefficienti_dir.size();
	long double C = C_rot_cin(coefficienti_dir);
	for (short i = 0 ; i < n  ; ++i){
		coefficienti_dir[i] /= C;
	}
	std::vector < long double > coefficienti;
	for (short i = 0 ; i < n ; ++i){
		coefficienti.push_back(coefficienti_dir[n - 1 - i]);
	}
		
	std::vector < long double > somma_geometrica;
	somma_geometrica.push_back(coefficienti[0]);
	for (short i = 1 ; i < n ; ++i){
		somma_geometrica.push_back(segno(somma_geometrica.back()) * sqrt(somma_geometrica.back()*somma_geometrica.back() + coefficienti[i]*coefficienti[i]));
	}
		
	std::vector< long double > fi;
	//NON RIPORTO IL 1° , è = 0 per definizione
	for (short i = 1 ; i < n  ; ++i){
		if((somma_geometrica[i] == 0.) && (coefficienti[i] == 0.)) {
			fi.push_back(0.);
		}
		else {
			long double fi_temp = coefficienti[i]/somma_geometrica[i];
			if (fi_temp > 1.) {
				if (fi_temp - 1. > errore_moltiplicazione_matrici) {
					std::cerr<<"ERRORE NEI COEFFICIENTI DI ROTAZIONE CINEMATICA" << std::endl;
				}
				else fi_temp = 1.;
			}
			else if (fi_temp < -1.) {
				if (-fi_temp - 1. > errore_moltiplicazione_matrici) {
					std::cerr<<"ERRORE NEI COEFFICIENTI DI ROTAZIONE CINEMATICA" << std::endl;
				}
				else fi_temp = -1.;
			}
			fi.push_back( std::acos(fi_temp));
		}
	}
	
	return fi;
	
}	

/*
 *	calcola gli angoli fi_ij del vettore di rotazione cinematica, dato il numero di particelle, i, j, nel caso di masse uguali, in modo k
 *
 */




inline std::vector<vettore_ldouble> rot_cin(const ublas::matrix<long double> & M_inversa) {
	unsigned short A = M_inversa.size1();
	std::vector< vettore_ldouble > rot_cin_temp;
	for (short i = 0; i < A ; ++i){
		for (short j = i+1; j < A; ++j){
			rot_cin_temp.push_back(rot_cin(M_inversa,i,j));
		}
	}
	unsigned short S = ((A-1)*A)/2;
	assert (rot_cin_temp.size() == S);
	return rot_cin_temp;
}


#endif
