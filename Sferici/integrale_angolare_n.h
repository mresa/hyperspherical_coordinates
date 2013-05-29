/*
 *  integrale_angolare_n.h
 *  Tesi
 *
 *  Created by Marco Resa on 21/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef INTEGRALE_ANGOLARE_N_H
#define INTEGRALE_ANGOLARE_N_H
#include "integrale_angolare.h"
#include "clebsch_n.h"

//integrale di armoniche sferiche per stati a n corpi , ridotto < l1,l2|L1><l3,L1|L2> …
inline long double int_nj_ridotto (const std::vector <unsigned short> & l, 
								   const std::vector <unsigned short> & L) {
	assert (l.size() == (L.size()+1));
	if (dis_triang (L, l)){
		long double coefficiente = int_3j_ridotto(l[0], l[1], L[0]);
		for (short i = 2 ; ( ( i < l.size() ) && ( coefficiente != 0.) ) ; ++i){		
			coefficiente *= int_3j_ridotto(l[i], L[i-2], L[i-1]);
		}
		return coefficiente;
	}
	else return 0.;
}

//integrale di armoniche sferiche per stati a n corpi , parte in m
inline long double int_nj_m (const std::vector <unsigned short> & l, 
							 const std::vector <short> & m, 
							 const std::vector <unsigned short> & L, 
							 const std::vector <short> & M) {
	controlla_l_m(l,m);
	controlla_l_m(L,M);
	assert (l.size() == (L.size()+1));
	long double coefficiente = int_3j_m(l[0], l[1], L[0], m[0], m[1], M[0]);
	for (short i = 2 ; ( ( i < l.size() ) && ( coefficiente != 0.) ) ; ++i){
		coefficiente *= int_3j_m(l[i], L[i-2], L[i-1], m[i], M[i-2], M[i-1]);
	}
	return coefficiente;
}

//integrale di armoniche sferiche per stati a n corpi 
inline long double int_nj (const std::vector <unsigned short> & l, 
						   const std::vector <short> & m, 
						   const std::vector <unsigned short> & L, 
						   const std::vector <short> & M) {
	controlla_l_m(l,m);
	controlla_l_m(L,M);
	assert (l.size() == (L.size()+1));
	long double risultato = int_nj_m(l,m,L,M);
	if (risultato != 0.) return int_nj_ridotto(l,L)*risultato;
	else return 0.;
}

//integrale di armoniche sferiche per stati a n corpi <[l],[m]|[l],[L],M>
inline long double int_nj (const std::vector <unsigned short> & l, 
						   const std::vector <short> & m, 
						   const std::vector <unsigned short> & L, 
						   const short & M) {
	assert (std::abs(M) <= L.back());
	short M_calcolato = 0.;
	for(short i = 0; i < m.size() ; ++i){		
		M_calcolato += (m[i]);
	}
	if (M_calcolato != M) return 0.;
	else return (int_nj(l, m, L, M_combinazioni_di_m (m)));
}


inline std::vector <long double> coefficienti_m (const std::vector <unsigned short> & l,
												 const std::vector <vettore_short> & m,
												 const std::vector <unsigned short> & L, 
												 const short & M) {
	std::vector <long double> coefficienti_m_temp;
	for(short i = 0; i < m.size() ; ++i){		
		coefficienti_m_temp.push_back(clebsch_gordan(l,m[i],L,M));
	}
	return coefficienti_m_temp;
}


inline std::vector <long double> coefficienti_m (const std::vector <unsigned short> & l,
												 const std::vector <unsigned short> & L, 
												 const short & M_A) {
	return coefficienti_m (l, costruisci_m(l , L, M_A) , L , M_A);
}

inline long double C_LMK (const std::vector <unsigned short> & l,
						  const std::vector <short> & m,
						  const unsigned short & L,
						  const short & M) {
	if(linear_sum(m) != M) return 0.;
	else{
		std::vector <vettore_ushort> BASE (somma_triangolare ( l, L ) );
		long double C = 0.;
		std::vector <short> M_generato (M_combinazioni_di_m (m));
#pragma omp parallel for reduction(+:C) schedule( guided )
		for(short i = 0; i < BASE.size() ; ++i){		
			long double C1 = int_nj_m( l, m, BASE[i], M_generato );
			if (C1 != 0) C1*= int_nj_ridotto(l, BASE[i]);
			C += C1;
		}
		return C;
	}
}



#endif

