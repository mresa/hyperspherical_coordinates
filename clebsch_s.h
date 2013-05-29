/*
 *  clebsch_s.h
 *  Tesi
 *
 *  Created by Marco Resa on 21/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef CLEBSCH_S_H
#define CLEBSCH_S_H
#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "algebra.h"
#include "triangolari_s.h"

inline long double clebsch_gordan_s (const unsigned short& dub_s1, 
									 const unsigned short& dub_s2, 
									 const unsigned short& dub_s,
									 const short& dub_z1, 
									 const short& dub_z2, 
									 const short& dub_z) {
	if (dub_z1 + dub_z2 != dub_z || dis_triang_s(dub_s,dub_s1,dub_s2) != true) return 0.;
	else {
		return (long double) (
							  pow((long double) -1., (long double)(dub_z+dub_s1-dub_s2)/2.)
							  *
							  sqrt((long double) dub_s + 1. )
							  *
							  gsl_sf_coupling_3j(dub_s1 , dub_s2 , dub_s ,dub_z1 ,dub_z2, - dub_z ));
	}
}

inline long double clebsch_gordan_s (const std::vector <unsigned short> & dub_s, 
									 const std::vector <short> & dub_z, 
									 const std::vector <unsigned short> & dub_S, 
									 const std::vector <short> & dub_Z) {
	assert (dub_s.size() == (dub_S.size()+1));
	long double coefficiente = 0.;
	if(controlla_s_z(dub_s,dub_z) && controlla_s_z(dub_S,dub_Z)){
		coefficiente = clebsch_gordan_s(dub_s[0], dub_s[1], dub_S[0], dub_z[0], dub_z[1], dub_Z[0]);
		
		for (short i = 2 ; ( ( i < dub_s.size() ) && ( coefficiente != 0.) ) ; ++i){
			coefficiente *= clebsch_gordan_s(dub_s[i], dub_S[i-2], dub_S[i-1], dub_z[i], dub_Z[i-2], dub_Z[i-1]);
		}
	}
	else;
	
	return coefficiente;
}

inline long double clebsch_gordan_s (const std::vector <unsigned short> & dub_s, 
									 const std::vector <short> & dub_z, 
									 const std::vector <unsigned short> & dub_S, 
									 const short & dub_Z) {
	assert (dub_s.size() == (dub_S.size()+1));
	long double coefficiente = 0.;
	if(controlla_s_z(dub_S.back(), dub_Z) && controlla_s_z(dub_s,dub_z)){
		short M_calcolato = 0.;
		for(short i = 0; i < dub_z.size() ; ++i){		
			M_calcolato += (dub_z[i]);
		}
		if (M_calcolato != dub_Z);
		else {
			std::vector <short> M_comb(M_combinazioni_di_m(dub_z));
			
			if (controlla_s_z(dub_S, M_comb)){
				coefficiente = (clebsch_gordan_s(dub_s, dub_z, dub_S, M_comb));
			}
		}
	}
	else;
	return coefficiente;
}

inline long double clebsch_gordan_s_cluster_2 (const std::vector <unsigned short> & dub_s, 
											   const std::vector <short> & dub_z, 
											   const std::vector <unsigned short> & dub_S1, 
											   const std::vector <unsigned short> & dub_S2,
											   const unsigned short& dub_S,
											   const short & dub_Z) {
	long double coefficiente = 0.;
	if(controlla_s_z(dub_S, dub_Z) && controlla_s_z(dub_s,dub_z)){
		short M_calcolato = 0.;
		for(short i = 0; i < dub_z.size() ; ++i){		
			M_calcolato += (dub_z[i]);
		}
		if (M_calcolato != dub_Z);
		else {
			std::vector <unsigned short> dub_s1;
			std::vector <unsigned short> dub_s2;
			std::vector <short> dub_z1;
			short Z1 = 0;
			
			std::vector <short> dub_z2;
			short Z2 = 0;
			for (short i = 0 ; i<=dub_S1.size(); ++i) {
				dub_s1.push_back(dub_s[i]);
				dub_z1.push_back(dub_z[i]);
				Z1 += dub_z1.back();
			}
			for (short i = 0 ; i<=dub_S2.size(); ++i) {
				dub_s2.push_back(dub_s[i+dub_z1.size()]);
				dub_z2.push_back(dub_z[i+dub_z1.size()]);
				Z2 += dub_z2.back();
			}
			
			if(controlla_s_z(dub_S1.back(), Z1) && controlla_s_z(dub_S2.back(), Z2)){
				coefficiente += clebsch_gordan_s(dub_S1.back(),dub_S2.back(),dub_S,Z1,Z2,dub_Z) * clebsch_gordan_s(dub_s1,dub_z1,dub_S1,Z1) *clebsch_gordan_s(dub_s2,dub_z2,dub_S2,Z2);
			}
		}
	}
	else;
	return coefficiente;
}

inline long double clebsch_gordan_s_cluster (const std::vector <vettore_ushort> & dub_s, 
											 const std::vector <vettore_short> & dub_z, 
											 const std::vector <vettore_ushort> & dub_S_r, 
											 const std::vector <unsigned short> & dub_S,
											 const short & dub_Z) {
	long double coefficiente = 0.;
	
	if(controlla_s_z(dub_S.back(), dub_Z) && controlla_s_z(dub_s,dub_z)){
		short M_calcolato = 0.;
		std::vector <short> dub_Z_r;
		std::vector <unsigned short> dub_S_r_ultimi;

		for(short spazio = 0; spazio < dub_z.size() ; ++spazio){		
			short dub_Z_r_t=0;
			for(short i = 0; i < dub_z[spazio].size() ; ++i){		
				M_calcolato += (dub_z[spazio][i]);
				dub_Z_r_t += (dub_z[spazio][i]);
			}
			dub_Z_r.push_back(dub_Z_r_t);
			dub_S_r_ultimi.push_back(dub_S_r[spazio].back());			
		}
		if (M_calcolato != dub_Z);
		else {			
			if(controlla_s_z(dub_S_r_ultimi, dub_Z_r) ){
				long double coefficiente_temp = 1;
				for(short spazio = 0; spazio < dub_z.size() ; ++spazio){		
					coefficiente_temp *= clebsch_gordan_s(dub_s[spazio],dub_z[spazio],dub_S_r[spazio],dub_Z_r[spazio]);
				}
				
				coefficiente += clebsch_gordan_s(dub_S_r_ultimi,dub_Z_r,dub_S,dub_Z) * coefficiente_temp;
				
			}
		}
	}
	else;
	return coefficiente;
}

#endif
