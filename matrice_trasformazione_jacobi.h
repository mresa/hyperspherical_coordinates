/*
 *  matrice_trasformazione_jacobi
 *  Tesi
 * 
 INPUT MASSE, MASSA_REF, cluster, 
 *
 *
 *	NOTA BENE : SCELTA PER L'ORDINAMENTO PROGRESSIVO DELLE COORDINATE DI JACOBI
 *	NEL CASO DI CLUSTER n_1, … , n_A si ottiene , con R_i centro massa del cluster i
 
 n_1 - 1 coordinate jacobi modo k
 1 coordinata relativa centro di massa 1 - centro di massa 2
 n_2 - 1 coordinate jacobi modo k
 1 coordinata relativa centro di massa (1,2) - centro di massa 3
 …
 …
 n_A - 1 coordinate jacobi modo k
 1 coordinata centro di massa totale (1,2,…,A)
 
 *  Created by Marco Resa on 24/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef MATRICE_JACOBI_H
#define MATRICE_JACOBI_H

#include <assert.h>
#include <vector>
#include "algebra.h"
#include "algebra_matrici.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
/*
 *	matrice di trasformazione per le coordinate di jacobi in forma vettoriale , per A particelle, nel modo K , con masse uguali
 *	M r = x in normalizzazione diretta
 */

inline ublas::matrix<long double> matrice_diretta(const unsigned short& A) {
	assert(A > 0);
	ublas::banded_matrix<long double> M(A,A,A,1);
	
	for(short i = 0 ; i < A - 1 ; ++i){
		//Fisso l'ultima riga
		M(A - 1,i) = 1./A; 
		//e la diagonale superiore
		M(i,i+1) = sqrt((long double) 2. * (i + 1.) / ( i + 2. ) );

		for(short j = 0 ; j <= i ; ++j)	
			M(i,j) = - sqrt((long double) 2. / ( (i + 1.) * ( i + 2. ) ) );
	}	
	//fisso l'ultimo elemento
	M(A - 1,A - 1) = 1./A;
	return M;
}

/*
 *	matrice di trasformazione per le coordinate di jacobi in forma vettoriale , per A particelle, nel modo K , con masse diverse
 *	M r = x
 */

inline ublas::matrix<long double> matrice_diretta(const std::vector < long double >& MASSE,  
												  const long double& MASSA_RIFERIMENTO) {
	assert(MASSE.size() >= 1);
	assert(MASSE[0] > 0.);
	assert(MASSA_RIFERIMENTO > 0.);
	
	unsigned short A = MASSE.size();
	
	std::vector < long double > M_j;
	M_j.push_back(MASSE[0]);
	
	for(short i = 1; i < A ; ++i) {
		assert( MASSE[i] > 0.);
		M_j.push_back(M_j.back() + MASSE[i]);
	}
		
	ublas::banded_matrix<long double> M(A,A,A,1);
	
	for(short i = 0 ; i < A - 1 ; ++i){
		
		M(i,i+1) = sqrt((long double) 2. * MASSE[i + 1.] * M_j[i] / ( MASSA_RIFERIMENTO * M_j[i + 1] )  );
		
		M( A - 1,i) = MASSE[i]/M_j.back();
		
		for(short j = 0 ; j <= i ; ++j){
			M(i,j) = - sqrt((long double) 2. * MASSE[i + 1.] * M_j[i] / 
							( MASSA_RIFERIMENTO * M_j[i + 1] ) ) * MASSE[j] / M_j[i];
		}
	}	
	M(A - 1, A - 1 ) = MASSE.back()/M_j.back();
	
	return M;
}

/*
 *	matrice di trasformazione per le coordinate di jacobi in forma vettoriale , per A particelle con masse diverse
 *	M r = x
 * CLUSTER: possibili modi H tramite input : es. A = 4 cluster = [2,2] -> modo H ; cluster = [1,1,1,1] -> modo k oppure cluster = [4] -> modo k
 */

inline ublas::matrix<long double> matrice_diretta(const std::vector < long double >& MASSE,  
												  const long double& MASSA_RIFERIMENTO, 
												  const std::vector < unsigned short >& CLUSTER){	
	int A = 0;
	
	for(short i = 0; i < CLUSTER.size() ; ++i) {
		assert (CLUSTER[i]>0.);
		A += CLUSTER[i];
	}
	assert (A = MASSE.size());
	
	
	MATRICI MATRICE(ublas::zero_matrix<long double>(A,A));
		
	//popolo la matrice temporanea con le matrici a blocchi
	
	unsigned short blocco = 0;
	
	for(short n = 0; n < CLUSTER.size() ; ++n) {
		
		std::vector < long double > MASSE_TEMP;
		
		unsigned short A_TEMP = CLUSTER[n];
		
		for(short m = blocco; m < blocco + A_TEMP ; ++m) {
			MASSE_TEMP.push_back(MASSE[m]);
		}
		
		MATRICI MATRICE_TEMP(matrice_diretta(MASSE_TEMP, MASSA_RIFERIMENTO));
				
		for(short i = 0; i < A_TEMP ; ++i) {			
			for(short j = 0; j < A_TEMP ; ++j) {
				MATRICE(blocco+i,blocco+j) = MATRICE_TEMP(i,j);
			}
		}
		
		blocco += A_TEMP;		
	}
		
	unsigned short A_1 = CLUSTER[0];
	long double MASSA_1 = 0.;
	for(short n = 0 ; n < A_1 ; ++n) {
		MASSA_1 += MASSE[n];
	}
	
	for(unsigned short n = 1 ; n < CLUSTER.size() ; ++n) {
		std::vector <long double> MASSE_TEMP;
		long double MASSA_2 = 0.;
		unsigned short A_2 = CLUSTER[n];
		for(short m = A_1 ; m < A_1 + A_2 ; ++m) {
			MASSA_2 += MASSE[m];
		}
		
		MASSE_TEMP.push_back(MASSA_1);
		MASSE_TEMP.push_back(MASSA_2);
		
		MATRICI MATRICE_TEMP(ublas::identity_matrix<long double>(A,A));
		MATRICI M_DIR(matrice_diretta(MASSE_TEMP,MASSA_RIFERIMENTO));
		
		MATRICE_TEMP(A_1 - 1,A_1 - 1) = M_DIR(0,0);
		MATRICE_TEMP(A_1 - 1,A_1 + A_2 - 1) = M_DIR(0,1);
		MATRICE_TEMP(A_1 + A_2 - 1,A_1 - 1) = M_DIR(1,0);
		MATRICE_TEMP(A_1 + A_2 - 1,A_1 + A_2 - 1) = M_DIR(1,1);
		A_1 += A_2;
		MASSA_1 += MASSA_2;
		MATRICE = ublas::prod_mkl(MATRICE_TEMP, MATRICE);
	}

	return MATRICE;
}

/*
 *	matrice di trasformazione per le coordinate di jacobi in forma vettoriale , per A particelle  con masse uguali alla massa di riferimento
 *	M r = x
 * CLUSTER: possibili modi H tramite input : es. A = 4 cluster = [2,2] -> modo H
 */

inline ublas::matrix<long double> matrice_diretta(const std::vector < unsigned short >& CLUSTER){	
	int A = 0;
	for(short i = 0; i < CLUSTER.size() ; ++i) {
		assert (CLUSTER[i]>0.);
		A += CLUSTER[i];
	}
	std::vector < long double > MASSE (A,1.0);	
	long double MASSA_RIFERIMENTO = 1.;
	return matrice_diretta(MASSE, MASSA_RIFERIMENTO, CLUSTER);
}


/*
 *	effetto sulla matrice della permutazione r_i <-> r_j
 *	NON SO PERCHÉ L'HO MESSA!!!!!!!!!!!
 *
 */

inline ublas::matrix<long double> Permutazione_ij(const ublas::matrix<long double> &M, 
												  const unsigned short& i, 
												  const unsigned short& j) {
	ublas::matrix<long double> MM(M);
	ublas::column (MM, j) = ublas::column (M, i);
	ublas::column (MM, i) = ublas::column (M, j);
	return MM;
}

#endif
