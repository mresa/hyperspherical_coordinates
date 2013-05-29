/*
 *  accoppiamenti1.h
 *  Tesi
 *
 *  Created by Marco Resa on 26/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ACC_ANG_H
#define ACC_ANG_H



//	i coefficienti di clebsh gordan per < j1,m1;j2,m2 |J2,M2>< j3,m3;J2,M2 |J3,M3> …… < jA,mA;J(A-1),M(A-1) |J_A,M_A> , in doppia precisione
// indico con L maiuscolo i relativi!!!!

#include <assert.h>
#include <iostream>
#include <iterator>
#include <vector>
#include <algorithm>
#include "algebra.h"
#include "combinazioni.h"
#include <limits>
#include <cmath>
#include "integrale_angolare_n.h"
#include "tipi_definiti.h"


/*
 *	restituisce lo sviluppo nella base [l] [m] dato un vettore della base [l][L],M
 *
 */

class Accoppiamenti {
	unsigned short NN;											// N totale
	std::vector < unsigned short > ll;							
	std::vector < unsigned short > LL;							
	short MM;													// M totale
	unsigned short degenerazione;								// numero di stati totali
	std::vector < long double > coefficiente;
	//std::vector < long double > coefficiente_int;
	std::vector < long double > coefficiente_lmk;
	std::vector < vettore_short > mm;							// vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	
	
public:
	
	
	Accoppiamenti (const std::vector < unsigned short >&  l , 
				   const std::vector < unsigned short >& L ,
				   const short&  M) {
		if (l.size() == 1){
			assert (l.front()==L.front()); //blocca lo script, da corregere?
			NN = l.size();
			ll = l;
			LL = ll;
			MM = M;
			degenerazione = 1;
			std::vector <short> m_temp(1,MM);
			mm.push_back(m_temp);
			coefficiente.push_back(1.);
			coefficiente_lmk.push_back(0.); //non utilizzato
											//coefficiente_int.push_back(0.); //non utilizzato
		}
		else{
			assert(l.size() > 1); // almeno 2 o non funziona
			assert(std::abs(M) <= L.back());
			assert(dis_triang(L, l)); //blocca lo script, da correggere
		
			NN = l.size();
			LL = L;
			MM = M;
			ll = l;
			degenerazione = 0;
		
			std::vector <vettore_short> M_temp (costruisci_M (LL , MM));
			
			for(short i = -ll[0] ; i <= ll[0] ; ++i){
				for(short j = 0 ; j < M_temp.size() ; ++j){
					const std::vector <short> m_temp (m_combinazioni_di_M(M_temp[j] , i));
					if (controlla_l_m(ll,m_temp)){
						long double coeff = clebsch_gordan (ll, m_temp, LL, MM);
						if (coeff != 0) {
							mm.push_back(m_temp);
							coefficiente.push_back(coeff);
							//coefficiente_int.push_back( int_nj_m_C(ll, m_temp, LL, M_temp[j]) * int_nj_ridotto_C(ll, LL) );
							coefficiente_lmk.push_back( C_LMK(ll, m_temp, LL.back() , MM));
						}
						else;
					}
					else;
				}
			}
			degenerazione = coefficiente.size();
		}
	}
	
	//distruttore
	~Accoppiamenti () {
		std::vector < vettore_short >().swap(mm);		
		std::vector < long double >().swap(coefficiente);
		std::vector < unsigned short >().swap(ll);
		std::vector < unsigned short >().swap(LL);
	}
	
	// N
	unsigned short N() const {return NN;}
	
	// L totale
	std::vector< unsigned short > L() const {return LL;}
	unsigned short L(const unsigned short& i) const {
		assert( i < NN - 1  );
		return LL[i];
	}	
	unsigned short L_tot() const {return LL.back();}
	
	// M totale
	short M() const {return MM;}
	
	
	// accesso agli l
	std::vector< unsigned short > l() const {return ll;}
	unsigned short l(const unsigned short& i) const {
		assert( i < NN );
		return ll[i];
	}	
	
	unsigned short D() const {return degenerazione;}	
	
	// accesso agli m
	std::vector< vettore_short > m() const {return mm;}
	std::vector< short > m(const unsigned short& i) const {
		assert( i < degenerazione );
		return mm[i];
	}
	short m(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < NN );
		return mm[i][j];
	}	
	
	// accesso ai coefficienti
	std::vector< long double > cm() const {
		return coefficiente;
	}
	long double cm(const unsigned short& i) const{
		assert( i < degenerazione );
		return coefficiente[i];
	}
	// accesso ai coefficienti
	/*std::vector< long double > ci() const {
		return coefficiente_int;
	}
	long double ci(const unsigned short& i) const{
		assert( i < degenerazione );
		return coefficiente_int[i];
	}
*/	// accesso ai coefficienti
	std::vector< long double > clmk() const {
		return coefficiente_lmk;
	}
	long double clmk(const unsigned short& i) const{
		assert( i < degenerazione );
		return coefficiente_lmk[i];
	}
};

#endif
