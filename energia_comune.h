/*
 *  energia_comune.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/01/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef ENERGIA_COMUNE_H
#define ENERGIA_COMUNE_H

std::vector < unsigned short > conversione_indici_isospin (const std::vector < short > & dub_zz, const std::vector < short > & dub_zz_ref){
	unsigned short A = dub_zz.size();
	std::vector < short > dub_zz_ref_t(dub_zz_ref);
	std::vector < unsigned short > temp;
	for (short i = 0 ; i<A; ++i) {
		temp.push_back(i);
	}
	for (short i=0; i< A; ++i) {
		if(dub_zz_ref_t[i]!=dub_zz[i]){ //se sono diversi lo cerco avanti e becco il primo che trovo
			for (short j = i ; j<A; ++j) {
				if(dub_zz_ref_t[j]==dub_zz[i]){ //se sono diversi lo cerco avanti
					dub_zz_ref_t[j] = dub_zz_ref[i];
					dub_zz_ref_t[i] = dub_zz[i];
					temp[i] = j;
					temp[j] = i;
				}
			}
		}
	}
	return temp;
}

inline void correzione_errore(MATRICI_SPARSE * M){
	for (short i = 0 ; i < (*M).size1() ; ++i){
		for (short j = i ; j < (*M).size2() ; ++j){
			if (std::abs((*M)(i,j)) <= err_matrici ) {
				(*M).erase_element(i,j);
				(*M).erase_element(j,i);
			}
			else;
		}
	}
}

inline void correzione_errore_spinta(MATRICI * M){
	for (short i = 0 ; i < (*M).size1() ; ++i){
		for (short j = 0 ; j < (*M).size2() ; ++j){
			if ((*M)(i,j) < -0.8) {
				(*M)(i,j) = -1.;
			}
			else if ((*M)(i,j) > 0.8) {
				(*M)(i,j) = 1.;
			}
			else {
				(*M).erase_element(i,j);
			}
		}
	}
}

inline void correzione_errore(V_MATRICI_SPARSE * M){
	for (short n = 0; n < (*M).size(); ++n) {
		for (short i = 0 ; i < (*M)[n].size1() ; ++i){
			for (short j = i ; j < (*M)[n].size2() ; ++j){
				if ( std::abs((*M)[n](i,j)) <= err_matrici ) {
					(*M)[n].erase_element(i,j);
					(*M)[n].erase_element(j,i);
				}
				else;
			}
		}
	}
}

inline void correzione_errore(MATRICI_SIMMETRICHE_UP* M){
	for (short i = 0 ; i < (*M).size1() ; ++i){
		for (short j = i ; j < (*M).size2() ; ++j){
			if ( std::abs((*M)(i,j)) <= err_matrici ) (*M)(i,j) = 0.;
			else;
		}
	}
}

inline void correzione_errore(MATRICI* M){
	for (short i = 0 ; i < (*M).size1() ; ++i){
		for (short j = 0 ; j < (*M).size2() ; ++j){
			if ( std::abs((*M)(i,j)) <= err_matrici) (*M)(i,j) = 0.;
			else;
		}
	}
}

inline long double correzione_errore(const long double & LD){
	if ( std::abs(LD) <= err_matrici ) return 0.;
	else return LD;
}

inline double correzione_errore(const double & D){
	if ( std::abs(D) <= err_matrici ) return 0.;
	else return D;
}
#endif
