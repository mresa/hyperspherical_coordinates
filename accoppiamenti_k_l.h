/*
 *  accoppiamenti_L.h
 *  Tesi
 *
 *  Created by Marco Resa on 11/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ACC_ANG_K_H
#define ACC_ANG_K_H

#include "accoppiamenti_angolari.h"

/*
 * 
 * dato L, M, K_l = l1 + … + ln genera gli stati (l1, … ln; L2, … , L; M) e ne fa lo sviluppo nella base l1,m1,l2,m2 … ln, mn con i relativi coefficienti
 * 
 * 
 */

class Accoppiamenti_k_l {
	
	unsigned short NN;											//
	unsigned short LL;											
	short MM;													// M totale
	unsigned short K_l;
	unsigned short degenerazione_base;							// numero di stati della base [l],[L],M
	unsigned short degenerazione_l;								// numero di stati della base [l],[m]    ?????
	std::vector <unsigned short> degenerazione_lm;				// numero di stati della base [l],[m] per ogni stato [l],[L],M
	std::vector < vettore_ushort > ll;							// vettore dei valori
	std::vector < vettore_ushort > LL_rel;						// numeri quantici L (N-1)
	std::vector < v_vettore_short > mm;							// vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	std::vector < vettore_ldouble > coefficiente;
	//std::vector < vettore_ldouble > coefficiente_int;
	std::vector < vettore_ldouble > coefficiente_lmk;


	
	//genero gli stati tali che l1 + … + ln = K_l
	static inline std::vector < vettore_ushort > init_l(const unsigned short& N, 
														const unsigned short& K_L) { //Ltotale
		std::vector <unsigned short> l_quantici(N, 0);
		std::vector <vettore_ushort> Arr_l;
		
		l_quantici.back() = K_L;			
		do {				
			Arr_l.push_back(l_quantici);
		} while (combinazioni(l_quantici.begin(), l_quantici.end()));
		
		return 	Arr_l;
	}
	//genero gli stati [L] tali che L_N = L  a partire dal vettore di stati l1 … ln e lo scrivo se L è giusto 
	static inline std::vector < v_vettore_ushort > init_L(const std::vector < vettore_ushort > & l,
														  const unsigned short & L) { //Ltotale
		std::vector <v_vettore_ushort> Arr_L;
		for(short i = 0; i < l.size(); ++i){
			std::vector < vettore_ushort > temp_somma (somma_triangolare(l[i]));
			std::vector < vettore_ushort > temp_L;
			for(short j = 0; j < temp_somma.size(); ++j){
				if (temp_somma[j].back() == L) temp_L.push_back(temp_somma[j]);
				else temp_L.push_back(temp_somma[j]);
			}
			
			if (temp_L.size() != 0) Arr_L.push_back(temp_L);
			else;
		}
		return Arr_L;
	}
	//riorganizzo la base
	static inline void init_base_l_L(const std::vector < vettore_ushort > & l,
									 const std::vector < v_vettore_ushort > & L,
									 std::vector < vettore_ushort > * l_base,
									 std::vector < vettore_ushort > * L_base) {
		for(short i = 0 ; i < l.size() ; ++i){
			for(short j = 0 ; j < L[i].size() ; ++j){
				(*l_base).push_back(l[i]);
				(*L_base).push_back(L[i][j]);
			}
		}
	}
		
	static inline void init_base(const unsigned short & N,
								 const unsigned short & L,
								 const unsigned short & K_L,
								 std::vector < vettore_ushort > * l_base,
								 std::vector < vettore_ushort > * L_base) {
		std::vector < vettore_ushort > l_temp(init_l (N,K_L));
		std::vector < v_vettore_ushort > L_temp (init_L(l_temp, L));
		for(short i = 0 ; i < l_temp.size() ; ++i){
			for(short j = 0 ; j < L_temp[i].size() ; ++j){
				if (L_temp[i][j].back() == L) {
					(*l_base).push_back(l_temp[i]);
					(*L_base).push_back(L_temp[i][j]);
				}
				else;
			}
		}
	}
	
	
	//devo costruire uno sviluppo per ogni vettore della base

public:
	Accoppiamenti_k_l ( const unsigned short  &  N, const unsigned short  &  L, const short&  M, const unsigned short&  K_L) {
		assert(std::labs((short)M) <= L);	
		LL = L;
		NN = N;
		K_l = K_L;
		MM = M;
		degenerazione_base = 0;
		degenerazione_l = 0;
		if (NN == 1) {
			if (K_l == L) {
				vettore_ushort l_temp;
				l_temp.push_back(L);
				ll.push_back(l_temp);
				LL_rel.push_back(l_temp);
				degenerazione_base = 1;
				Accoppiamenti TEMP(l_temp, l_temp, MM);
				if (TEMP.D()!=0) {
					degenerazione_l += TEMP.D();
					degenerazione_lm.push_back(TEMP.D());
					mm.push_back(TEMP.m());
					coefficiente.push_back(TEMP.cm());
					coefficiente_lmk.push_back(TEMP.clmk()); //non utilizzato
				}
				else;
				
			}
		}
		else if (NN > 1) {
			init_base( NN, LL, K_l, &ll, &LL_rel);
			degenerazione_base = ll.size();
			for(short i = 0; i < degenerazione_base; ++i){
				Accoppiamenti TEMP(ll[i], LL_rel[i], MM);
				if (TEMP.D()!=0) {
					degenerazione_l += TEMP.D();
					degenerazione_lm.push_back(TEMP.D());
					mm.push_back(TEMP.m());
					coefficiente.push_back(TEMP.cm());
					//coefficiente_int.push_back(TEMP.ci());
					coefficiente_lmk.push_back(TEMP.clmk());
				}
				else;
			}			
		}
	}	
	
	
	//distruttore
	~Accoppiamenti_k_l () {
	}
	
	inline unsigned short N() const {return NN;}
	inline unsigned short L_TOT() const {return LL;}											
	inline short M() const {return MM;}
	inline unsigned short K_L() const {return K_l;}
	inline unsigned short D() const {return degenerazione_base;}	
	inline unsigned short D_TOT() const {return degenerazione_l;}	

	inline unsigned short DLM(const unsigned short& i) const {
		assert( i < degenerazione_base );
		return degenerazione_lm[i];
	}
	
	inline std::vector< vettore_ushort > l() const {return ll;}
	inline std::vector< unsigned short > l(const unsigned short& i) const {
		assert( i < degenerazione_base );
		return ll[i];
	}
	
	inline std::vector< vettore_ushort > L() const {return LL_rel;}
	inline std::vector< unsigned short > L(const unsigned short& i) const {
		assert( i < degenerazione_base );
		return LL_rel[i];
	}
	
	
	inline std::vector< v_vettore_short > m() const {return mm;}
	inline std::vector< vettore_short > m(const unsigned short& i) const {
		assert( i < degenerazione_base );
		return mm[i];
	}
	inline std::vector< short > m(const unsigned short& i, const unsigned short& j) const {
		assert( i < degenerazione_base );
		assert( j < degenerazione_lm[i] );
		return mm[i][j];
	}
	inline short m(const unsigned short& i, const unsigned short& j, const unsigned short& k) const {
		assert( i < degenerazione_base );
		assert( j < degenerazione_lm[i] );
		assert( k < NN );
		return mm[i][j][k];
	}
	
	inline std::vector <vettore_ldouble> cm() const {
		return coefficiente;
	}

	inline std::vector< long double > cm(const unsigned short& i) const {
		assert( i < degenerazione_base );		
		return coefficiente[i];
	}
	
	inline long double cm(const unsigned short& i, const unsigned short& j) const{
		assert( i < degenerazione_base );
		assert( j < degenerazione_lm[i] );
		return coefficiente[i][j];
	}
	
/*	inline std::vector <vettore_ldouble> ci() const {
		return coefficiente_int;
	}
	
	inline std::vector< long double > ci(const unsigned short& i) const {
		assert( i < degenerazione_base );		
		return coefficiente_int[i];
	}
	
	inline long double ci(const unsigned short& i, const unsigned short& j) const{
		assert( i < degenerazione_base );
		assert( j < degenerazione_lm[i] );
		return coefficiente_int[i][j];
	}
*/	
	inline std::vector <vettore_ldouble> clmk() const {
		return coefficiente_lmk;
	}
	
	inline std::vector< long double > clmk(const unsigned short& i) const {
		assert( i < degenerazione_base );		
		return coefficiente_lmk[i];
	}
	
	inline long double clmk(const unsigned short& i, const unsigned short& j) const{
		assert( i < degenerazione_base );
		assert( j < degenerazione_lm[i] );
		return coefficiente_lmk[i][j];
	}
	
};

#endif
