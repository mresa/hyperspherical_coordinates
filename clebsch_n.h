/*
 *  clebsh_n.h
 *  Tesi
 *
 *  Created by Marco Resa on 21/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef CLEBSCH_N_H
#define CLEBSCH_N_H

#include "clebsch.h"
/*
 *	calcolo dei coefficienti a molti corpi
 */


// < [l],[L],[M]|[l],[m]>

inline long double clebsch_gordan (const std::vector <unsigned short> & l, 
								   const std::vector <short> & m, 
								   const std::vector <unsigned short> & L, 
								   const std::vector <short> & M){
	if( controlla_l_m(l,m) && controlla_l_m(L,M)){
		assert (l.size() == (L.size()+1));
		long double coefficiente = clebsch_gordan(l[0], l[1], L[0], m[0], m[1], M[0]);
		for (short i = 2 ; ( ( i < l.size() ) && ( coefficiente != 0.) ) ; ++i){
			coefficiente *= clebsch_gordan(l[i], L[i-2], L[i-1], m[i], M[i-2], M[i-1]);
		}
		return coefficiente;
	}
	else return 0;
}

// < [l],[L],M|[l],[m]>
inline long double clebsch_gordan (const std::vector <unsigned short> & l, 
								   const std::vector <short> & m, 
								   const std::vector <unsigned short> & L, 
								   const short & M) {
	assert (std::abs(M) <= L.back());
	assert (l.size() == (L.size()+1));
	short M_calcolato = 0.;
	for(short i = 0; i < m.size() ; ++i){		
		M_calcolato += (m[i]);
	}
	if (M_calcolato != M) return 0.;
	else return (clebsch_gordan(l, m, L, M_combinazioni_di_m (m)));
}

#endif

