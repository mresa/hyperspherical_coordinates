/*
 *  accoppiamenti_iperangolari.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/06/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ACC_IPER_H
#define ACC_IPER_H

#include "accoppiamenti_k_l.h"
#include "stringhe.h"

/*
 * 
 * dato N , L , M e K
 *
 * genera la base [n2 … nN], [l1, … ,lN], [L2, … L], M
 * ne fa lo sviluppo nella base [n2 … nN], [l1, … ,lN] ,[m1 … mN] a K fisso
 * 
 * 
 */

class Accoppiamenti_Iperangolari {
	
	unsigned short NN;												// N
	unsigned short KK;												// K totale
	unsigned short LL;												// L totale
	short MM;														// M totale
	unsigned long degenerazione;									// degenerazione della base
	std::vector <unsigned long> degenerazione_m;					// degenerazione della parte m, per ogni stato con [n_i][l_i] fisso, indicizzato da i
	std::vector < vettore_ushort > nn;								// vettore dei valori [n] : NN[i] = [n]
	std::vector < vettore_ushort > ll;								// vettore dei valori [l] : LL[i] = [l]
	std::vector < vettore_ushort > LL_rel;							// vettore dei valori [l] : LL[i] = [l]
	std::vector < v_vettore_short > mm;								// vettore dei valori [m] : MM[i][j] = [m] relativi allo stato NN[i] LL[i]
	std::vector < vettore_ldouble > coefficiente;					// vettore dei coefficienti per lo sviluppo nella base
	//std::vector < vettore_ldouble > coefficiente_int;					// vettore dei coefficienti per lo sviluppo nella base
	std::vector < vettore_ldouble > coefficiente_lmk;					// vettore dei coefficienti per lo sviluppo nella base
	std::vector < vettore_bool > parita_stato;						// vettore dei coefficienti per lo sviluppo nella base

	
	// ritorna un vettore con le possibili combinazioni tali che somma n_i = K_n
	static std::vector < vettore_ushort > init_N (const unsigned short& N, 
												  const unsigned short& K_N ) {
		std::vector <unsigned short> n_quantici(N-1, 0);
		std::vector <vettore_ushort> ArrN;
		
		n_quantici.back() = K_N;			
		do {				
			ArrN.push_back(n_quantici);				
		} while (combinazioni(n_quantici.begin(), n_quantici.end()));
		
		return 	ArrN;
	}
	
	static inline bool controlla_parita(const unsigned short& P, const std::vector < unsigned short > & CLUSTER, const std::vector < unsigned short > &l){
		unsigned short pp = (P %2);
		unsigned short N = l.size();
		
		bool verita = (l.back() %2) == pp;
		
		for (short i = 0 ; i < CLUSTER.size() - 1 ; ++i) {
			verita *= (l [  N - 1 - CLUSTER [ CLUSTER.size() - 1  - i ] ] %2 == pp);
		}
		return verita;
	}
	
	static std::vector < bool > genera_parita_stato (const std::vector < unsigned short > & CLUSTER, const std::vector < unsigned short > &l){
		std::vector < bool > parita;
		unsigned short N = l.size();
		parita.push_back(l.back()%2);
		for (short i = 0 ; i < CLUSTER.size() - 1 ; ++i) {
			parita.push_back(l [  N - 1 - CLUSTER [ CLUSTER.size() - 1  - i ] ] %2 );
		}
		return parita;
	}
	
	static std::vector < bool > genera_parita_stato ( const std::vector < unsigned short > &l){
		std::vector < bool > parita;
		parita.push_back(l.back()%2);
		return parita;
	}

public:
	
	Accoppiamenti_Iperangolari (const unsigned short & N, 
								const unsigned short & L,
								const short & M, 
								const unsigned short & K) {
		assert(N > 0);
		assert(std::labs(M) <= L);
		NN = N;
		KK = K;
		LL = L;
		MM = M;
		degenerazione = 0.;
		if (NN == 1) {
			Accoppiamenti_k_l TEMP( NN, LL, MM, KK);
			std::vector <unsigned short> n_quantici;
			for(short j = 0; j < TEMP.D() ; ++j){
				nn.push_back(n_quantici); // ho un vettore vuoto!
				ll.push_back(TEMP.l(j));
				parita_stato.push_back(genera_parita_stato(ll.back()));
				LL_rel.push_back(TEMP.L(j));
				mm.push_back(TEMP.m(j));
				coefficiente.push_back(TEMP.cm(j));
				coefficiente_lmk.push_back(TEMP.clmk(j));
				degenerazione_m.push_back(TEMP.DLM(j));
				++degenerazione;
			}
		}
		else {
			for (short kl = KK ; kl >= 0; kl -= 2 ) {
				std::vector <vettore_ushort> ArrN( init_N (NN, (KK-kl)/2) );
				Accoppiamenti_k_l TEMP( NN, LL, MM, kl);
				
				for (short i = 0; i < ArrN.size() ; ++i) {
					for(short j = 0; j < TEMP.D() ; ++j){
						nn.push_back(ArrN[i]);
						ll.push_back(TEMP.l(j));
						parita_stato.push_back(genera_parita_stato(ll.back()));
						LL_rel.push_back(TEMP.L(j));
						mm.push_back(TEMP.m(j));
						coefficiente.push_back(TEMP.cm(j));
						//coefficiente_int.push_back(TEMP.ci(j));
						coefficiente_lmk.push_back(TEMP.clmk(j));
						degenerazione_m.push_back(TEMP.DLM(j));
						++degenerazione;
					}
				}
			}
		}
	}
	
	Accoppiamenti_Iperangolari (const unsigned short & N, 
								const unsigned short & L,
								const short & M,
								const unsigned short & K, 
								const unsigned short& P, 
								const std::vector < unsigned short > & CLUSTER ) {
		assert(N > 0);
		assert(std::labs(M) <= L);
		NN = N;
		KK = K;
		LL = L;
		MM = M;
		degenerazione = 0.;
		if(P==2){
			if (NN == 1) {
				Accoppiamenti_k_l TEMP( NN, LL, MM, KK);
				std::vector <unsigned short> n_quantici;
				for(short j = 0; j < TEMP.D() ; ++j){
					nn.push_back(n_quantici); // ho un vettore vuoto!
					ll.push_back(TEMP.l(j));
					parita_stato.push_back(genera_parita_stato(ll.back()));
					LL_rel.push_back(TEMP.L(j));
					mm.push_back(TEMP.m(j));
					coefficiente.push_back(TEMP.cm(j));
					coefficiente_lmk.push_back(TEMP.clmk(j));
					degenerazione_m.push_back(TEMP.DLM(j));
					++degenerazione;
				}
			}
			else{
				for (short kl = KK ; kl >= 0; kl -= 2 ) {
					std::vector <vettore_ushort> ArrN( init_N (NN, (KK-kl)/2) );
					Accoppiamenti_k_l TEMP( NN, LL, MM, kl);
					for (short i = 0; i < ArrN.size() ; ++i) {
						for(short j = 0; j < TEMP.D() ; ++j){
							
							nn.push_back(ArrN[i]);
							ll.push_back(TEMP.l(j));
							parita_stato.push_back(genera_parita_stato(CLUSTER,ll.back()));
							LL_rel.push_back(TEMP.L(j));
							mm.push_back(TEMP.m(j));
							coefficiente.push_back(TEMP.cm(j));
							coefficiente_lmk.push_back(TEMP.clmk(j));
							degenerazione_m.push_back(TEMP.DLM(j));
							++degenerazione;
						}
					}
				}
			}
		}
		else{
			if (NN == 1) {
				Accoppiamenti_k_l TEMP( NN, LL, MM, KK);
				std::vector <unsigned short> n_quantici;
				for(short j = 0; j < TEMP.D() ; ++j){
					if ( controlla_parita( P,CLUSTER,TEMP.l(j) ) ) {

					nn.push_back(n_quantici); // ho un vettore vuoto!
					ll.push_back(TEMP.l(j));
					parita_stato.push_back(genera_parita_stato(ll.back()));
					LL_rel.push_back(TEMP.L(j));
					mm.push_back(TEMP.m(j));
					coefficiente.push_back(TEMP.cm(j));
					coefficiente_lmk.push_back(TEMP.clmk(j));
					degenerazione_m.push_back(TEMP.DLM(j));
					++degenerazione;
					}
				}
			}
			else{
				for (short kl = KK ; kl >= 0; kl -= 2 ) {
					std::vector <vettore_ushort> ArrN( init_N (NN, (KK-kl)/2) );
					Accoppiamenti_k_l TEMP( NN, LL, MM, kl);
					
					for (short i = 0; i < ArrN.size() ; ++i) {
						for(short j = 0; j < TEMP.D() ; ++j){
							if ( controlla_parita( P,CLUSTER,TEMP.l(j) ) ) {
								nn.push_back(ArrN[i]);
								ll.push_back(TEMP.l(j));
								parita_stato.push_back(genera_parita_stato(CLUSTER,ll.back()));
								LL_rel.push_back(TEMP.L(j));
								mm.push_back(TEMP.m(j));
								coefficiente.push_back(TEMP.cm(j));
								//coefficiente_int.push_back(TEMP.ci(j));
								coefficiente_lmk.push_back(TEMP.clmk(j));
								degenerazione_m.push_back(TEMP.DLM(j));
								++degenerazione;
							}
							else;
						}
					}
				}
			}
		}
	}
	
	//distruttore
	~Accoppiamenti_Iperangolari () {
	}
	
	// N° di particelle -1
	unsigned short N() const {return NN;}
	
	unsigned short L_TOT() const {return LL;}											

	short M() const {return MM;}

	// K totale
	unsigned short K() const {return KK;}
	
	// degenerazione
	unsigned long D() const {return degenerazione;}
	std::vector < unsigned long > DM() const {return degenerazione_m;}
	unsigned long DM(const unsigned long& i) const {
		assert( i < degenerazione );
		return degenerazione_m[i];
	}
	
	
	// vettore dei valori [n], per lo stato i-esimo
	
	std::vector < vettore_ushort > n() const {return nn;}
	std::vector< unsigned short > n(const unsigned long& i) const {
		assert( i < degenerazione );
		return nn[i];
	}
	
	// vettore dei valori [l], per lo stato i-esimo
	std::vector < vettore_ushort > l() const {return ll;}
	std::vector< unsigned short > l(const unsigned long& i) const {
		assert( i < degenerazione );
		return ll[i];
	}
	std::vector < vettore_bool > P_l() const {return parita_stato;}
	std::vector< bool > P_l(const unsigned long& i) const {
		assert( i < degenerazione );
		return parita_stato[i];
	}
	
	
	// vettore dei valori [L], per lo stato i-esimo
	std::vector < vettore_ushort > L() const {return LL_rel;}
	std::vector< unsigned short > L(const unsigned long& i) const {
		assert( i < degenerazione );
		return LL_rel[i];
	}
	
	// vettore dei valori [m], per lo stato i-esimo in n,l; stato j-esimo della degenerazione
	std::vector< v_vettore_short > m() const {return mm;}
	std::vector< vettore_short > m(const unsigned long& i) const {
		assert( i < degenerazione );
		return mm[i];
	}
	std::vector< short > m(const unsigned long& i, const unsigned long& j) const {
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return mm[i][j];
	}
	
	// vettore dei coefficienti per lo stato [m] per avere L,M fisso
	std::vector <vettore_ldouble> cm() const {
		return coefficiente;
	}	
	std::vector< long double > cm(const unsigned long& i) const {
		assert( i < degenerazione );		
		return coefficiente[i];
	}	
	long double cm(const unsigned long& i, const unsigned long& j) const{
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return coefficiente[i][j];
	}		


/*	std::vector <vettore_ldouble> ci() const {
		return coefficiente_int;
	}	
	std::vector< long double > ci(const unsigned long& i) const {
		assert( i < degenerazione );		
		return coefficiente_int[i];
	}	
	long double ci(const unsigned long& i, const unsigned long& j) const{
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return coefficiente_int[i][j];
	}		

*/	std::vector <vettore_ldouble> clmk() const {
		return coefficiente_lmk;
	}	
	std::vector< long double > clmk(const unsigned long& i) const {
		assert( i < degenerazione );		
		return coefficiente_lmk[i];
	}	
	long double clmk(const unsigned long& i, const unsigned long& j) const{
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return coefficiente_lmk[i][j];
	}		
	
};

/*
 * 
 * dato L , M e K_max ne fa lo sviluppo nella base n2, l1,m1,l2,m2, fino a K_max
 * 
 * 
 */

class Accoppiamenti_Iperangolari_K_max {
		
	unsigned short NN;												
	unsigned short KK_max;											// K totale
	unsigned short LL;												// L totale
	short MM;														// M totale
	unsigned long degenerazione;									// degenerazione della base
	std::vector <unsigned long> degenerazione_m;					// degenerazione della parte m, per ogni stato con [n_i][l_i] fisso, indicizzato da i
	std::vector < vettore_ushort > nn;								// vettore dei valori [n] : NN[i] = [n]
	std::vector < vettore_ushort > ll;								// vettore dei valori [l] : LL[i] = [l]
	std::vector < vettore_ushort > LL_rel;								
	std::vector < v_vettore_short > mm;								// vettore dei valori [m] : MM[i][j] = [m] relativi allo stato NN[i] LL[i]
	std::vector < unsigned short > KK;								// valore di K per ogni stato
	std::vector < vettore_ldouble > coefficiente;					// vettore dei coefficienti per lo sviluppo nella base
	//std::vector < vettore_ldouble > coefficiente_int;					// vettore dei coefficienti per lo sviluppo nella base
	std::vector < vettore_bool > parita_stato;					// vettore dei coefficienti per lo sviluppo nella base
	unsigned short parita_spaziale;
	std::string s;
	std::string sd;
public:
	

	Accoppiamenti_Iperangolari_K_max (const unsigned short&  N ,
									  const unsigned short &  L,
									  const short&  M,
									  const unsigned short&  K_M,
									  const unsigned short& P,
									  const std::vector < unsigned short > & CLUSTER, 
									  const unsigned short& parita ) {
		assert(std::labs(M) <= L);
		NN = N;
		KK_max = K_M;
		LL = L;
		MM = M;
		degenerazione = 0.;
		parita_spaziale=parita;
		unsigned short partenza;
		unsigned short incremento = 2;
		if ( parita_spaziale == 0 ) partenza = 0;		//parità totale 1
		else if( parita_spaziale == 1 ) partenza = 1;	//parità totale -1
		else if( parita_spaziale == 2 ) {				//tutte le parità
			partenza = 0;
			incremento = 1;
		}
		for(unsigned short k = partenza; k <= KK_max; k+=incremento){
			Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k,P,CLUSTER);
			
			for (short i = 0; i < TEMP.D() ; ++i) {
				nn.push_back(TEMP.n(i));
				ll.push_back(TEMP.l(i));
				parita_stato.push_back(TEMP.P_l(i));
				LL_rel.push_back(TEMP.L(i));
				coefficiente.push_back(TEMP.cm(i));
				//coefficiente_int.push_back(TEMP.ci(i));
				mm.push_back(TEMP.m(i));
				degenerazione_m.push_back(TEMP.DM(i));
				++degenerazione;
				KK.push_back(k);
			}
		}		
		s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(P) + "_Par" + scrivi(parita_spaziale) + "]";
		sd = "N" + scrivi(NN) + "__P_" + scrivi(P) + "__Par_" + scrivi(parita_spaziale);
	}
	
	
	Accoppiamenti_Iperangolari_K_max (const unsigned short&  N ,
									  const unsigned short &  L,
									  const short&  M,
									  const unsigned short&  K_M,
									  const unsigned short& P,
									  const std::vector < unsigned short > & CLUSTER ) {
		assert(std::labs(M) <= L);
		NN = N;
		KK_max = K_M;
		LL = L;
		MM = M;
		degenerazione = 0.;
		parita_spaziale = 2;
		for(short k = 0; k <= KK_max; k+=2){
			Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k,P,CLUSTER);
			
			for (short i = 0; i < TEMP.D() ; ++i) {
				nn.push_back(TEMP.n(i));
				ll.push_back(TEMP.l(i));
				parita_stato.push_back(TEMP.P_l(i));
				LL_rel.push_back(TEMP.L(i));
				coefficiente.push_back(TEMP.cm(i));
				//coefficiente_int.push_back(TEMP.ci(i));
				mm.push_back(TEMP.m(i));
				degenerazione_m.push_back(TEMP.DM(i));
				++degenerazione;
				KK.push_back(k);
			}
		}	
		s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(P) + "_Par" + scrivi(parita_spaziale) + "]";
		sd = "N" + scrivi(NN) + "__P_" + scrivi(P) + "__Par_" + scrivi(parita_spaziale);

	}
	
	Accoppiamenti_Iperangolari_K_max (const unsigned short&  N ,
									  const unsigned short &  L, 
									  const short&  M,
									  const unsigned short&  K_M ,
									  const unsigned short& parita) {
		assert(std::labs(M) <= L);
		NN = N;
		KK_max = K_M;
		LL = L;
		MM = M;
		degenerazione = 0.;
		parita_spaziale = parita;
		unsigned short partenza;
		unsigned short incremento = 2;
		if ( parita_spaziale == 0 ) partenza = 0;		//parità totale 1
		else if( parita_spaziale == 1 ) partenza = 1;	//parità totale -1
		else if( parita_spaziale == 2 ) {				//tutte le parità
			partenza = 0;
			incremento = 1;
		}
		for(short k = partenza; k <= KK_max; k+=incremento){
			Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k);
			
			for (short i = 0; i < TEMP.D() ; ++i) {
				nn.push_back(TEMP.n(i));
				ll.push_back(TEMP.l(i));
				parita_stato.push_back(TEMP.P_l(i));
				LL_rel.push_back(TEMP.L(i));
				coefficiente.push_back(TEMP.cm(i));
				//coefficiente_int.push_back(TEMP.ci(i));
				mm.push_back(TEMP.m(i));
				degenerazione_m.push_back(TEMP.DM(i));
				++degenerazione;
				KK.push_back(k);
			}
			
		}	
		s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(2) + "_Par" + scrivi(parita_spaziale) + "]";
		sd = "N" + scrivi(NN) + "__P_" + scrivi(2) + "__Par_" + scrivi(parita_spaziale);

	}
	
	Accoppiamenti_Iperangolari_K_max (const unsigned short&  N ,
									  const unsigned short &  L, 
									  const short&  M,
									  const unsigned short&  K_M ) {
		assert(std::labs(M) <= L);
		NN = N;
		KK_max = K_M;
		LL = L;
		MM = M;
		degenerazione = 0.;
		parita_spaziale = 2;
		for(short k = 0; k <= KK_max; k+=2){ 
			// FIXME: suppongo pari ?? avevo usato k+=2 forse cos' non tornano le dimensioni della base
			// TODO: K pari poiché la parità di k me la da la parità di L, che userò pari
			Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k);
			
			for (short i = 0; i < TEMP.D() ; ++i) {
				nn.push_back(TEMP.n(i));
				ll.push_back(TEMP.l(i));
				parita_stato.push_back(TEMP.P_l(i));
				LL_rel.push_back(TEMP.L(i));
				coefficiente.push_back(TEMP.cm(i));
				//coefficiente_int.push_back(TEMP.ci(i));
				mm.push_back(TEMP.m(i));
				degenerazione_m.push_back(TEMP.DM(i));
				++degenerazione;
				KK.push_back(k);
				
			}
		}	
		s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(2) + "_Par" + scrivi(parita_spaziale) + "]";
		sd = "N" + scrivi(NN) + "__P_" + scrivi(2) + "__Par_" + scrivi(parita_spaziale);
	}
	
	
	//distruttore
	~Accoppiamenti_Iperangolari_K_max () {
	}
	
	//aggiungono stati!!
/*	void PLUS(const unsigned short&  K_plus){
		if(K_plus > 0){
			for(short k = KK_max + 1; k <= KK_max + K_plus ; ++k){ 
				Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k);
				
				for (short i = 0; i < TEMP.D() ; ++i) {
					nn.push_back(TEMP.n(i));
					ll.push_back(TEMP.l(i));
					parita_stato.push_back(TEMP.P_l(i));
					LL_rel.push_back(TEMP.L(i));
					coefficiente.push_back(TEMP.cm(i));
					mm.push_back(TEMP.m(i));
					degenerazione_m.push_back(TEMP.DM(i));
					++degenerazione;
					KK.push_back(k);
					
				}
			}
			KK_max += K_plus;
		}
	}
	
	void PLUS(const unsigned short&  K_plus,
			  const unsigned short& P,
			  const std::vector < unsigned short > & CLUSTER){
		if(K_plus > 0){
			for(short k = KK_max + 1; k <= KK_max + K_plus; ++k){
				Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k,P,CLUSTER);
				
				for (short i = 0; i < TEMP.D() ; ++i) {
					nn.push_back(TEMP.n(i));
					ll.push_back(TEMP.l(i));
					parita_stato.push_back(TEMP.P_l(i));
					LL_rel.push_back(TEMP.L(i));
					coefficiente.push_back(TEMP.cm(i));
					mm.push_back(TEMP.m(i));
					degenerazione_m.push_back(TEMP.DM(i));
					++degenerazione;
					KK.push_back(k);
				}
			}
			KK_max += K_plus;
		}
	}
*/	
	void PLUS(const unsigned short& K_plus){
		if(K_plus > 0){
			unsigned short partenza;
			unsigned short incremento = 2;
			if ( parita_spaziale == 0 && KK_max %2 == 0) partenza = KK_max + 2;		//parità totale 1
			else if ( parita_spaziale == 0 && KK_max %2 == 1) partenza = KK_max + 1;		//parità totale 1
			else if( parita_spaziale == 1 && KK_max %2 == 0) partenza = KK_max + 1;						//parità totale -1
			else if( parita_spaziale == 1 && KK_max %2 == 1) partenza = KK_max + 2;						//parità totale -1
			else if( parita_spaziale == 2 ) {									//tutte le parità
				partenza = KK_max + 1;
				incremento = 1;
			}
			for(short k = partenza; k <= KK_max + K_plus; k+=incremento){
				Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k);
				
				for (short i = 0; i < TEMP.D() ; ++i) {
					nn.push_back(TEMP.n(i));
					ll.push_back(TEMP.l(i));
					parita_stato.push_back(TEMP.P_l(i));
					LL_rel.push_back(TEMP.L(i));
					coefficiente.push_back(TEMP.cm(i));
					//coefficiente_int.push_back(TEMP.ci(i));
					mm.push_back(TEMP.m(i));
					degenerazione_m.push_back(TEMP.DM(i));
					++degenerazione;
					KK.push_back(k);
				}
			}
			KK_max += K_plus;
			s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(2) + "_Par" + scrivi(parita_spaziale) + "]";
			sd = "N" + scrivi(NN) + "__P_" + scrivi(2) + "__Par_" + scrivi(parita_spaziale);
		}
	}
	
	void PLUS(const unsigned short&  K_plus, 
			  const unsigned short& P, 
			  const std::vector < unsigned short > & CLUSTER){
		if(K_plus > 0){
			unsigned short partenza;
			unsigned short incremento = 2;
			if ( parita_spaziale == 0 && KK_max %2 == 0) partenza = KK_max + 2;		//parità totale 1
			else if ( parita_spaziale == 0 && KK_max %2 == 1) partenza = KK_max + 1;		//parità totale 1
			else if( parita_spaziale == 1 && KK_max %2 == 0) partenza = KK_max + 1;						//parità totale -1
			else if( parita_spaziale == 1 && KK_max %2 == 1) partenza = KK_max + 2;						//parità totale -1
			else if( parita_spaziale == 2 ) {									//tutte le parità
				partenza = KK_max + 1;
				incremento = 1;
			}
			for(unsigned short k = partenza; k <= KK_max + K_plus; k+=incremento){
				Accoppiamenti_Iperangolari TEMP(NN,LL,MM,k,P,CLUSTER);
				for (unsigned long i = 0; i < TEMP.D() ; ++i) {
					nn.push_back(TEMP.n(i));
					ll.push_back(TEMP.l(i));
					parita_stato.push_back(TEMP.P_l(i));
					LL_rel.push_back(TEMP.L(i));
					coefficiente.push_back(TEMP.cm(i));
					mm.push_back(TEMP.m(i));
					degenerazione_m.push_back(TEMP.DM(i));
					++degenerazione;
					KK.push_back(k);
				}
			}
			KK_max += K_plus;
			s = "A" + scrivi(NN + 1) + "_K_min_" + scrivi(KK_max) + "_[P" + scrivi(P) + "_Par" + scrivi(parita_spaziale) + "]";
			sd = "N" + scrivi(NN) + "__P_" + scrivi(P) + "__Par_" + scrivi(parita_spaziale);
		}
	}
	
	
	inline std::string filename() const { return s; }
	
	inline std::string fileD() const { 
		std::string sd1 = "_D_" + scrivi(degenerazione);
		std::string sdo = sd + sd1;
		return sdo; 
	}
	inline std::string fileD(const unsigned long & d) const { 
		std::string sd1 = "_D_" + scrivi(d);
		std::string sdo = sd + sd1;
		return sdo; 
	}

	
	// N° di particelle -1
	inline unsigned short N() const {return NN;}
	
	inline unsigned short P() const {return parita_spaziale;}
	
 	inline unsigned short L_TOT() const {return LL;}											
	
	inline short M() const {return MM;}
	
	// K totale
	inline unsigned short K_M() const {return KK_max;}
	
	// degenerazione
	inline unsigned long D() const {return degenerazione;}
	inline std::vector < unsigned long > DM() const {return degenerazione_m;}
	inline unsigned long DM(const unsigned long& i) const {
		assert( i < degenerazione );
		return degenerazione_m[i];
	}
	
	
	// vettore dei valori [n], per lo stato i-esimo
	
	inline std::vector < vettore_ushort > n() const {return nn;}
	inline std::vector< unsigned short > n(const unsigned long& i) const {
		assert( i < degenerazione );
		return nn[i];
	}
	
	inline std::vector < unsigned short > k() const {return KK;}
	inline unsigned short k(const unsigned long& i) const {
		assert( i < degenerazione );
		return KK[i];
	}

	// vettore dei valori [l], per lo stato i-esimo
	inline std::vector < vettore_ushort > l() const {return ll;}
	inline std::vector< unsigned short > l(const unsigned long& i) const {
		assert( i < degenerazione );
		return ll[i];
	}
	
	std::vector < vettore_bool > P_l() const {return parita_stato;}
	std::vector< bool > P_l(const unsigned long& i) const {
		assert( i < degenerazione );
		return parita_stato[i];
	}
	
	// vettore dei valori [L], per lo stato i-esimo
	inline std::vector < vettore_ushort > L() const {return LL_rel;}
	inline std::vector< unsigned short > L(const unsigned long& i) const {
		assert( i < degenerazione );
		return LL_rel[i];
	}
	
	// vettore dei valori [m], per lo stato i-esimo in n,l; stato j-esimo della degenerazione
	inline std::vector< v_vettore_short > m() const {return mm;}
	inline std::vector< vettore_short > m(const unsigned long& i) const {
		assert( i < degenerazione );
		return mm[i];
	}
	inline std::vector< short > m(const unsigned long& i, const unsigned long& j) const {
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return mm[i][j];
	}
	
	// vettore dei coefficienti per lo stato [m] per avere L,M fisso
	inline std::vector <vettore_ldouble> cm() const {
		return coefficiente;
	}	
	inline std::vector< long double > cm(const unsigned long& i) const {
		assert( i < degenerazione );		
		return coefficiente[i];
	}	
	inline long double cm(const unsigned long& i, const unsigned long& j) const{
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return coefficiente[i][j];
	}



/*	// vettore dei coefficienti per lo stato [m] per avere L,M fisso
	inline std::vector <vettore_ldouble> ci() const {
		return coefficiente_int;
	}	
	inline std::vector< long double > ci(const unsigned long& i) const {
		assert( i < degenerazione );		
		return coefficiente_int[i];
	}	
	inline long double ci(const unsigned long& i, const unsigned long& j) const{
		assert( i < degenerazione );
		assert( j < degenerazione_m[i] );
		return coefficiente_int[i][j];
	}
*/
};





inline MATRICI_SPARSE L_12(const Accoppiamenti_Iperangolari_K_max & TEMP){
	MATRICI_SPARSE MATRI(TEMP.D(),TEMP.D());
	for (short i = 0 ; i < TEMP.D() ; ++i){
		long double L = TEMP.l(i).back();
		MATRI(i,i) = L;
	}
	return MATRI;
}
inline MATRICI_SPARSE L_12_parita(const Accoppiamenti_Iperangolari_K_max & TEMP,
								  const unsigned short P){
	MATRICI_SPARSE MATRI(TEMP.D(),TEMP.D());
	for (short i = 0 ; i < TEMP.D() ; ++i){
		unsigned short L = TEMP.l(i).back();
		if (L%2==P){
			MATRI(i,i) = 1;
		}
	}
	return MATRI;
}


#endif