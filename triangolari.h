/*
 *  triangolari.h
 *  Tesi
 *
 *  Created by Marco Resa on 25/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRIANGOLARI_H
#define TRIANGOLARI_H

#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "tipi_definiti.h"
#include "algebra.h"
#include "algebra_vettori.h"
#include "combinazioni.h"


/*
 *	Verifica che il vettore degli l e degli m siano compatibili (stessa dimensione e |m| < l)
 */

inline bool controlla_l_m (const unsigned short & l, 
						   const short & m){
	bool verita = (std::abs(m) <= l);
	return verita;
}

/*inline bool controlla_l_m (const std::vector <unsigned short> & l, 
						   const std::vector <short> & m){
	bool verita = (bool) (l.size() == m.size());
	for( short i = 0; verita * (i < l.size()) ; ++i){
		verita *= (std::abs(m[i]) <= l[i]);
	}
	return verita;
}
*/
inline bool controlla_l_m (const std::vector <unsigned short> & l, 
						   const std::vector <short> & m){
	bool verita = (bool) (l.size() == m.size());
#pragma omp parallel for reduction(*:verita) schedule( guided )
	for( short i = 0; i < l.size() ; ++i){
		verita *= (std::abs(m[i]) <= l[i]);
		if (verita == false) i = l.size();
	}
	return verita;
}

/*
 *	Disuguaglianza triangolare (iterativa)
 */

// l1 + l2 -> L
inline bool dis_triang (const unsigned short & L, const unsigned short & l1, const unsigned short & l2) {
	bool verita = false;
	if (std::labs( l1 - l2 ) <= L && L <= l1 + l2){ // check sulla parità: pari pari va in pari pari dispari in dispari
		if( (L%2) == ((l1 + l2)%2) ) verita = true;
		else;
	}
	else;	
	return verita;
}

// [l] -> L
inline bool dis_triang (const unsigned short &L, const std::vector <unsigned short>& ll) {
	assert(ll.size() > 0);
	if(ll.size() == 1) return ( L == ll[0]);
	else if(ll.size() == 2) return dis_triang(L, ll[0], ll[1]);
	else {
		long L_min = ll[0];
		long L_max = ll[0];
		for (short i = 1; i < ll.size() ; ++i) {
			if( ll[i] <= L_max && ll[i] >= L_min ) L_min = 0;
			else if( ll[i] > L_max ) L_min = ll[i] - L_max;
			else if( ll[i] < L_min ) L_min = L_min - ll[i];
			L_max += ll[i];
		}
		
		return (L_min <= L && L <= L_max);
	}
}

// [l] -> [L]
inline bool dis_triang (const std::vector <unsigned short> &L, const std::vector <unsigned short>& ll) {
	assert(ll.size() == L.size() + 1 );
	bool verita = dis_triang(L[0], ll[0], ll[1]);
	for (short i = 1 ; i< L.size() ; ++i ) {
		verita *= dis_triang(L[i], L[i - 1] , ll[i + 1] );
	}
	return verita;
}

// l1[i] + l2[i] -> L[i]
inline bool dis_triang (const std::vector <unsigned short> &L, const std::vector <unsigned short>& l1, const std::vector <unsigned short>& l2) {
	unsigned short n = L.size();
	assert(l1.size() == n );
	assert(l2.size() == n );
	bool verita = true;
	for (short i = 0 ; i < n && verita != false ; ++i ) {
		verita *= dis_triang(L[i], l1[i], l2[i] );
	}
	return verita;
}

// -m1[i] + m2[i] -> M[i]
inline bool dis_triang_m (const std::vector <short> & M, const std::vector <short>& m1, const std::vector <short>& m2) {
	unsigned short n = M.size();
	assert(m1.size() == n );
	assert(m2.size() == n );
	bool verita = true;
	for (short i = 0 ; i < n && verita != false ; ++i ) {
		verita *= ((-m1[i] + m2[i]) == M[i]);
	}
	return verita;
}



/*
 * gli do l1 e l2 mi da la sequenza di possibili L come vettore = (|l1-l2| ,… , l1 + l2)
 * 
 *	valido per 1 + 2 = (12)
 *	
 */

inline std::vector <unsigned short> somma_triangolare_2 (const unsigned short & l1, const unsigned short & l2){
	std::vector <unsigned short> L;
	for (short i = std::abs(l2 - l1); i <= l1 + l2; ++i){
		L.push_back(i);
	}
	return L;
}

/*
 * gli do una sequenza di L come vettore (|l1-l2| ,… , l1 + l2),  aggiunge l e mi da tutte le possibili L : (L, L')^1 …  (L, L')^N
 *	valido per (12) + 3 = [(12), (123)]
 */

inline std::vector <vettore_ushort> somma_triangolare_3 (const std::vector <unsigned short> & L, const unsigned short & l){
	std::vector <vettore_ushort> L_tot;
	for (short i = 0; i < L.size(); ++i){
		std::vector <unsigned short> L_sum (somma_triangolare_2 ( L[i] , l));
		for (short j = 0; j < L_sum.size(); ++j){
			std::vector <unsigned short> L_temp;
			L_temp.push_back(L[i]);
			L_temp.push_back(L_sum[j]);
			L_tot.push_back(L_temp);
		}
	}
	return L_tot;
}

/*
 * gli do l1 e l2 mi da la sequenza di possibili L come VETTORI DI VETTORE a un valore	(|l1-l2|)
 (	…	)
 (l1 + l2)
 */

inline std::vector <vettore_ushort> somma_triangolare (const unsigned short & l1, const unsigned short & l2){
	std::vector <vettore_ushort> L;
	for (short i = std::abs(l2 - l1); i <= l1 + l2; ++i){
		std::vector <unsigned short> L_temp;
		L_temp.push_back(i);
		L.push_back(L_temp);
	}
	return L;
}

/*
 * gli do una sequenza di L possibili come vettore ([L]^1, … , [L]^N),  aggiunge l e mi da tutte le possibili L : ([L], L')^1 …  ([L], L')^N'   
 */

inline std::vector <vettore_ushort> somma_triangolare (const std::vector <vettore_ushort> & L, const unsigned short & l){
	std::vector <vettore_ushort> L_tot;
	for (short i = 0; i < L.size(); ++i){
		std::vector <unsigned short> L_sum (somma_triangolare_2 ( L[i].back() , l));
		for (short j = 0; j < L_sum.size(); ++j){
			std::vector <unsigned short> L_temp (L[i]);
			L_temp.push_back(L_sum[j]);
			L_tot.push_back(L_temp);
		}
	}
	return L_tot;
}

/*
 * gli do una sequenza di tanti l e lui  trova tutti i possibili [L12; L123; … ;L]
 */

inline std::vector <vettore_ushort> somma_triangolare (const std::vector <unsigned short> & l){
	
	//condizione debole
	assert (l.size() > 0);
	
	if (l.size()==1) {
		std::cerr<< "attenzione : 2 corpi" << std::endl;
		std::vector <vettore_ushort> L;
		L.push_back(l);
		return L;
	} 
	
	//fine condizione debole
	
	else{
		std::vector <vettore_ushort> L(somma_triangolare(l[0], l[1]));	
		for (short i = 2; i < l.size(); ++i){
			L = somma_triangolare(L, l[i]);
		}
		return L;
	}
}

/*
 * gli do una sequenza di tanti l e lui  trova L12; L123; … ;L; restituisce solo con L fisso
 */

inline std::vector <vettore_ushort> somma_triangolare (const std::vector <unsigned short> & l, const unsigned short & L) {
	
	std::vector <vettore_ushort> L_temp (somma_triangolare (l));
	std::vector <vettore_ushort> L_tot;
	for (short i = 0 ; i < L_temp.size() ; ++i) {
		if (L_temp[i].back()==L) {
			L_tot.push_back(L_temp[i]);
		}
		else;
	}
	return L_tot;	
}


/*
 *	dati gli m, calcolo [M]  : M_0 = m0+m1 M_1 = M_0+m2 etc
 */

inline std::vector <short> M_combinazioni_di_m (const std::vector <short> & m) {
	assert (m.size() > 0); //almeno 2 ?
	std::vector <short> M;
	if (m.size() == 1) {
		M.push_back(m[0]);
	}
	else {
		M.push_back(m[0]+m[1]);
		for (short i = 2 ; i < m.size() ; ++i){
			M.push_back( M.back() + m[i]);
		}	
	}
	return M;
}

/*
 *	dati gli M e un m1, calcolo gli m
 */

inline std::vector <short> m_combinazioni_di_M (const std::vector <short> & M, const short & m1) {
	assert (M.size() > 0); //almeno 2 !!
	std::vector <short> m;
	m.push_back(m1);
	m.push_back( M[0] - m1);
	
	for ( short i = 1 ; i < M.size() ; ++i){
		m.push_back( M[i] - M[i-1]);
	}
	return m;
}


/*
 *	somma diretta di spazi
 */


inline std::vector <vettore_short> somma_diretta_M (const std::vector <vettore_short> & M_1, const short & M_2){
	std::vector <vettore_short> M;
	for (short i = 0; i< M_1.size(); ++i){
		std::vector <short>  M_temp(M_1[i]);
		M_temp.push_back(M_2);
		M.push_back(M_temp);
	}
	return M;
}

inline std::vector <vettore_short> somma_diretta_M (const std::vector <vettore_short> & M_1, const std::vector <short> & M_2){
	std::vector <vettore_short> M;
	for (short i = 0; i< M_1.size(); ++i){
		std::vector <short>  M_temp(M_1[i]);
		for (short j1 = 0; j1< M_2.size(); ++j1){
			M_temp.push_back(M_2[j1]);
		}
		M.push_back(M_temp);
	}
	return M;
}

inline std::vector <vettore_short> somma_diretta_M (const std::vector <vettore_short> & M_1, const std::vector <vettore_short> & M_2){
	std::vector <vettore_short> M;
	for (short i = 0; i< M_1.size(); ++i){
		for (short j = 0; j< M_2.size(); ++i){
			std::vector <short>  M_temp(M_1[i]);
			for (short j1 = 0; j1< M_2[j].size(); ++j1){
				M_temp.push_back(M_2[j][j1]);
			}
			M.push_back(M_temp);
		}
	}
	return M;
}

inline std::vector <v_vettore_short> somma_diretta_M (const std::vector <v_vettore_short> & M_1, const std::vector <vettore_short> & M_2){
	std::vector <v_vettore_short> M;
	for (short i = 0; i< M_1.size(); ++i){
		for (short j = 0; j< M_2.size(); ++i){
			std::vector <vettore_short> M_t(M_1[i]);
			M_t.push_back(M_2[j]);
			M.push_back(M_t);
		}
	}
	return M;
}

inline std::vector <v_vettore_ushort> somma_diretta_M (const std::vector <v_vettore_ushort> & M_1, const std::vector <vettore_ushort> & M_2){
	std::vector <v_vettore_ushort> M;
	for (short i = 0; i< M_1.size(); ++i){
		for (short j = 0; j< M_2.size(); ++j){
			std::vector <vettore_ushort> M_t(M_1[i]);
			M_t.push_back(M_2[j]);
			M.push_back(M_t);
		}
	}
	return M;
}

//costruisce tutti i possibili valori di [m] a partire dai valori di [l]

inline std::vector <vettore_short> costruisci_M (const std::vector <unsigned short> & L) {
	assert (L.size() > 0);
	
	std::vector <short> mm;
	
	for(short i = 0; i < L.size() ; ++i){		
		mm.push_back(L[i]);
	}
	
	std::vector <vettore_short> M;
	std::vector <short> mm_l(mm);	
	do{
		controlla_l_m(L , mm);
		M.push_back(mm);
	}while (genera_m (mm.begin(), mm.end(), mm_l.end()));
	
	return M;
}

//costruisce tutti i possibili valori di [m] a partire dai valori di [l], con ultimo m fissato
inline std::vector <vettore_short> costruisci_M (const std::vector <unsigned short> & L , const short & M_A) {
	assert (std::abs(M_A) <= L.back());
	
	std::vector <unsigned short> l(L);
	std::vector <vettore_short> M;
	if (l.size() == 1) {
		vettore_short m;
		m.push_back(M_A);
		M.push_back(m);
	}
	else{
		l.pop_back(); 
		M = somma_diretta_M(costruisci_M(l),M_A);
	}	
	return M;
}

//costruisce tutti i possibili valori di [m] a partire dai valori di [l] e compatibilmente a m_1 + m_2 = M_1   M_1 + m_3 = M_2 etc
inline std::vector <vettore_short> costruisci_m (const std::vector <unsigned short> & l , const std::vector <vettore_short> & M) {
	
	assert (M[0].size() == l.size() - 1);
	
	std::vector <vettore_short> m;
	
	for(short i = -l[0] ; i <= l[0] ; ++i){
		for(short j = 0 ; j < M.size() ; ++j){
			vettore_short m_temp (m_combinazioni_di_M(M[j] , i));
			if (controlla_l_m(l,m_temp)) m.push_back(m_temp);
			else;
		}
	}
	return m;
}

//costruisce tutti i possibili valori di [m] a partire dai valori di [l], [L] , M_totale
inline std::vector <vettore_short> costruisci_m (const std::vector <unsigned short> & l , const std::vector <unsigned short> & L, const short & M_A) {
	return costruisci_m( l , costruisci_M( L , M_A) );
}


#endif

