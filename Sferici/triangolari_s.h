/*
 *  triangolari_s.h
 *  Tesi
 
 
		PER LA PARTE DI SPIN LAVORO SU numeri quantici * 2 !!!
 
 *
 *  Created by Marco Resa on 21/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TRIANGOLARI_S_H
#define TRIANGOLARI_S_H

#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "algebra.h"
#include "algebra_vettori.h"
#include "combinazioni.h"


inline bool controlla_s_z (const unsigned short & dub_s, 
						   const short & dub_z){
	bool verita = (std::abs(dub_z) <= dub_s) && ( (dub_s - std::abs(dub_z) ) % 2 == 0);
	return verita;
}

inline bool controlla_s_z (const std::vector <unsigned short> & dub_s, 
						   const std::vector <short> & dub_z){
	bool verita = (bool) (dub_s.size() == dub_z.size());
	for( short i = 0; verita*(i < dub_s.size()) ; ++i){
		verita *= controlla_s_z(dub_s[i], dub_z[i]);
	}
	return verita;
}

inline bool controlla_s_z (const std::vector <vettore_ushort> & dub_s, 
						   const std::vector <vettore_short> & dub_z){
	bool verita = (bool) (dub_s.size() == dub_z.size());
	for( short spazio = 0; verita*(spazio < dub_s.size()) ; ++spazio){
		for( short i = 0; verita*(i < dub_s[spazio].size()) ; ++i){
		verita *= controlla_s_z(dub_s[spazio][i], dub_z[spazio][i]);
		}
	}
	return verita;
}

/*
 *	Disuguaglianza triangolare (iterativa)
 */

// s1 + s2 -> S
inline bool dis_triang_s (const unsigned short & dub_S, const unsigned short & dub_s1, const unsigned short & dub_s2) {
	bool verita = false;
	if (std::labs( dub_s1 - dub_s2 ) <= dub_S && dub_S <= dub_s1 + dub_s2){ // check sulla parità, lavoro in dub_l : pari pari va in pari pari dispari in dispari
		if( (dub_S%2) == ((dub_s1 + dub_s2)%2) ) verita = true;
		else;
	}
	else;	
	return verita;
}

// [s] -> S_tot
inline bool dis_triang_s (const unsigned short & dub_S, const std::vector <unsigned short>& dub_s) {
	assert(dub_s.size() > 0);
	if(dub_s.size() == 1) return ( dub_S == dub_s[0]);
	else if(dub_s.size() == 2) return dis_triang_s(dub_S, dub_s[0], dub_s[1]);
	else {
		long S_min = dub_s[0];
		long S_max = dub_s[0];
		for (short i = 1; i < dub_s.size() ; ++i) {
			if( dub_s[i] <= S_max && dub_s[i] >= S_min ) S_min = 0;
			else if( dub_s[i] > S_max ) S_min = dub_s[i] - S_max;
			else if( dub_s[i] < S_min ) S_min = S_min - dub_s[i];
			S_max += dub_s[i];
		}
		bool P = (S_min%2);
		assert ( P ==(S_max%2));
		return (S_min <= dub_S && dub_S <= S_max && (dub_S%2) == P);
	}
}

// [s] -> [S_12, S_123 , … , S_tot]
inline bool dis_triang_s (const std::vector <unsigned short> & dub_S, 
						  const std::vector <unsigned short> & dub_s) {
	assert(dub_s.size() == dub_S.size() + 1 );
	bool verita = dis_triang_s(dub_S[0], dub_s[0], dub_s[1]);
	for (short i = 1 ; verita *(i< dub_S.size()) ; ++i ) {
		verita *= dis_triang_s(dub_S[i], dub_S[i - 1] , dub_s[i + 1] );
	}
	return verita;
}

// [s] -> [S_12, S_123 , … , S_tot] su tanti spazi
inline bool dis_triang_s (const std::vector <vettore_ushort> & dub_S, 
						  const std::vector <vettore_ushort> & dub_s) {
	assert(dub_s.size() == dub_S.size());
	bool verita = true;
	for (short i = 0 ; verita *(i< dub_S.size()) ; ++i ) {
		verita *= dis_triang_s(dub_S[i], dub_s[i]);
	}
	return verita;
}
// [s] -> [S_12, S_123 , … , S_tot] su tanti spazi
inline bool dis_triang_s (const std::vector < unsigned short >& dub_S,
						  const std::vector < vettore_ushort >& dub_S_r) {
	assert(dub_S.size() == dub_S_r.size() - 1 );
	std::vector < unsigned short > dub_S_r_corretto;
	for (short i = 0 ;  i< dub_S_r.size() ; ++i ) {
		dub_S_r_corretto.push_back(dub_S_r[i].back());
	}
	return dis_triang_s(dub_S,dub_S_r_corretto);
}

// [s1] + [s2] -> [S] in modo indipendente
inline bool dis_triang_s (const std::vector <unsigned short> & dub_S, 
						  const std::vector <unsigned short>& dub_s1, 
						  const std::vector <unsigned short>& dub_s2) {
	unsigned short n = dub_S.size();
	assert(dub_s1.size() == n );
	assert(dub_s2.size() == n );
	bool verita = true;
	for (short i = 0 ; i < n && verita != false ; ++i ) {
		verita *= dis_triang_s(dub_S[i], dub_s1[i], dub_s2[i] );
	}
	return verita;
}

// -z1[i] + z2[i] -> Z[i]
inline bool dis_triang_z (const std::vector <short> & Z, const std::vector <short>& z1, const std::vector <short>& z2) {
	unsigned short n = Z.size();
	assert(z1.size() == n );
	assert(z2.size() == n );
	bool verita = true;
	for (short i = 0 ; i < n && verita != false ; ++i ) {
		verita *= ((-z1[i] + z2[i]) == Z[i]);
	}
	return verita;
}

//costruisce i valori z a partire dai valori di s
inline std::vector <vettore_short> costruisci_Z (const std::vector <unsigned short> & dub_S) {
	assert (dub_S.size() > 0);
	
	std::vector <short> zz;
	
	for(short i = 0; i < dub_S.size() ; ++i){		
		zz.push_back(dub_S[i]);
	}
	
	std::vector <vettore_short> Z;
	std::vector <short> zz_t(zz);	
	do{
		if (controlla_s_z(dub_S , zz)) Z.push_back(zz);
	}while (genera_m_dub (zz.begin(), zz.end(), zz_t.end()));
	
	return Z;
}

//costruisce i valori m a partire dai valori di l, con ultimo m fissato
inline std::vector <vettore_short> costruisci_Z (const std::vector <unsigned short> & dub_S , 
												 const short & dub_Z) {
	std::vector <vettore_short> Z;
	if( controlla_s_z(dub_S.back(),dub_Z)){
		
		std::vector <unsigned short> dub_S_t(dub_S);
		
		if (dub_S_t.size() == 1) {
			vettore_short z;
			z.push_back(dub_Z);
			Z.push_back(z);
		}
		else{
			dub_S_t.pop_back(); 
			Z = somma_diretta_M(costruisci_Z(dub_S_t),dub_Z);
		}	
	}
	return Z;
}

//costruisce i valori m a partire dai valori di l, con somma fissata
inline std::vector <vettore_short> costruisci_Z_tot (const std::vector <unsigned short> & dub_S , 
													 const short & dub_Z_tot) {
	std::vector <vettore_short> Z_t(costruisci_Z(dub_S));
	std::vector <vettore_short> Z;
	
	for (short i = 0 ; i < Z_t.size(); ++i) {
		if( linear_sum(Z_t[i]) == dub_Z_tot){
			Z.push_back(Z_t[i]);
		}
	}
	return Z;
}

//CLUSTER


inline std::vector <vettore_ushort> clusterizza (const std::vector <unsigned short> & dub_S, const std::vector < unsigned short > & CLUSTER) {
	std::vector <vettore_ushort> risultato;
	unsigned short j_temp = 0;
	for (short i = 0 ; i<CLUSTER.size(); ++i) {
		vettore_ushort dub_s_cl_t;
		for (short j = 0 ; j<CLUSTER[i]; ++j) {
			dub_s_cl_t.push_back(dub_S[j+j_temp]);
		}
		risultato.push_back(dub_s_cl_t);
		j_temp += dub_s_cl_t.size();
	}
	return risultato;
}
	
inline std::vector <vettore_short> clusterizza (const std::vector <short> & dub_S, const std::vector < unsigned short > & CLUSTER) {
	std::vector <vettore_short> risultato;
	unsigned short j_temp = 0;
	for (short i = 0 ; i<CLUSTER.size(); ++i) {
		vettore_short dub_s_cl_t;
		for (short j = 0 ; j<CLUSTER[i]; ++j) {
			dub_s_cl_t.push_back(dub_S[j+j_temp]);
		}
		risultato.push_back(dub_s_cl_t);
		j_temp += dub_s_cl_t.size();
	}
	return risultato;
}

inline std::vector <unsigned short> declusterizza (const std::vector <vettore_ushort> & dub_S) {
	std::vector < unsigned short > dub_ss_temp(dub_S[0]);
	for(short spazi = 1 ; spazi < dub_S.size() ; ++spazi){
		for(short indici = 0 ; indici < dub_S[spazi].size() ; ++indici){
			dub_ss_temp.push_back(dub_S[spazi][indici]);
		}
	}
	return dub_ss_temp;
}

inline std::vector <short> declusterizza (const std::vector <vettore_short> & dub_S) {
	std::vector <short > dub_ss_temp(dub_S[0]);
	for(short spazi = 1 ; spazi < dub_S.size() ; ++spazi){
		for(short indici = 0 ; indici < dub_S[spazi].size() ; ++indici){
			dub_ss_temp.push_back(dub_S[spazi][indici]);
		}
	}
	return dub_ss_temp;
}
//costruisce i valori z a partire dai valori di s
inline std::vector <v_vettore_short> costruisci_Z (const std::vector <vettore_ushort> & dub_S) {
	assert (dub_S.size() > 0);
	vettore_ushort dub_S_linea(dub_S[0]);

	vettore_ushort CLUSTER;
	CLUSTER.push_back(dub_S[0].size());

	for (short spazio = 1; spazio < dub_S.size(); ++spazio) {
		for (short indice = 0; indice < dub_S[spazio].size(); ++indice) {
			dub_S_linea.push_back(dub_S[spazio][indice]);
		}

		CLUSTER.push_back(dub_S[spazio].size());
	}

	std::vector <vettore_short> zz(costruisci_Z (dub_S_linea));

	std::vector <vettore_ushort> ss_clusterizzato(clusterizza(dub_S_linea,CLUSTER));

	std::vector <v_vettore_short> ZZ;
	for (short casi = 0; casi < zz.size(); ++casi) {
		std::vector <vettore_short> zz_clusterizzato(clusterizza(zz[casi],CLUSTER));
		ZZ.push_back(zz_clusterizzato);
	}
	return ZZ;
}

//costruisce i valori z a partire dai valori di S con z totale fissato
inline std::vector <v_vettore_short> costruisci_Z (const std::vector <vettore_ushort> & dub_S, 
												   const short & dub_Z) {
	std::vector <v_vettore_short> z_temp(costruisci_Z(dub_S));
	std::vector <v_vettore_short> z_ok;
	for (short i = 0 ; i < z_temp.size(); ++i) {
		short Z=0;
		for (short j = 0 ; j < z_temp[i].size(); ++j) {
			for (short k = 0 ; k < z_temp[i][j].size(); ++k) {
				Z+=z_temp[i][j][k];
			}
		}
		if (Z==dub_Z){
			z_ok.push_back(z_temp[i]);
		}
	}
	return z_ok;
}

//costruisce i valori m a partire dai valori di l, con ultimo m fissato
inline std::vector <v_vettore_short> costruisci_Z_b (const std::vector <vettore_ushort> & dub_S, 
												   const short & dub_Z) {
	std::vector <v_vettore_short> zz (costruisci_Z (dub_S));
	std::vector <v_vettore_short> ZZ;
	for(short i = 0; i < zz.size() ; ++i){
		short somma = 0;
		for(short spazio = 0; spazio < zz[i].size() ; ++spazio){
			for(short indice = 0; indice < zz[i][spazio].size() ; ++indice){
				somma += zz[i][spazio][indice];
			}
		}
		if( somma == dub_Z ){
			ZZ.push_back(zz[i]);
		}
	}
	
	return ZZ;
}


/// FINE CLUSTER



/*
 * gli do l1 e l2 mi da la sequenza di possibili L come vettore = (|l1-l2| ,… , l1 + l2)
 * 
 *	valido per 1 + 2 = (12)
 *	
 */
inline std::vector <unsigned short> somma_triangolare_2_s (const unsigned short & dub_s1, 
														   const unsigned short & dub_s2){
	std::vector <unsigned short> dub_S;

	for (short i = std::abs(dub_s2 - dub_s1); i <= dub_s2 + dub_s1 ; i+=2){
		dub_S.push_back(i);
	}
	return dub_S;
}
/*
 * gli do l1 e l2 mi da la sequenza di possibili L come VETTORI DI VETTORE a un valore	(|l1-l2|)
 (	…	)
 (l1 + l2)
 */

inline std::vector <vettore_ushort> somma_triangolare_s (const unsigned short & dub_s1, 
														 const unsigned short & dub_s2){
	std::vector <vettore_ushort> dub_S;
	for (short i = std::abs(dub_s2 - dub_s1); i <= dub_s2 + dub_s1; i+=2){
		std::vector <unsigned short> dub_S_temp;
		dub_S_temp.push_back(i);
		dub_S.push_back(dub_S_temp);
	}
	return dub_S;
}

/*
 * gli do una sequenza di L possibili come vettore ([L]^1, … , [L]^N),  aggiunge l e mi da tutte le possibili L : ([L], L')^1 …  ([L], L')^N'   
 */

inline std::vector <vettore_ushort> somma_triangolare_s (const std::vector <vettore_ushort> & dub_S, 
														 const unsigned short & dub_s){
	std::vector <vettore_ushort> dub_S_tot;
	for (short i = 0; i < dub_S.size(); ++i){
		std::vector <unsigned short> dub_S_sum (somma_triangolare_2_s ( dub_S[i].back(), dub_s));
		for (short j = 0; j < dub_S_sum.size(); ++j){
			std::vector <unsigned short> dub_S_temp (dub_S[i]);
			dub_S_temp.push_back(dub_S_sum[j]);
			dub_S_tot.push_back(dub_S_temp);
		}
	}
	return dub_S_tot;
}

/*
 * gli do una sequenza di tanti l e lui  trova tutti i possibili [L12; L123; … ;L]
 */

inline std::vector <vettore_ushort> somma_triangolare_s (const std::vector <unsigned short> & dub_s){
	
	//condizione debole
	assert (dub_s.size() > 0);
	
	if (dub_s.size()==1) {
		std::cerr<< "attenzione : 2 corpi" << std::endl;
		std::vector <vettore_ushort> dub_S;
		dub_S.push_back(dub_s);
		return dub_S;
	} 
	
	//fine condizione debole
	
	else{
		std::vector <vettore_ushort> dub_S(somma_triangolare_s(dub_s[0], dub_s[1]));	
		for (short i = 2; i < dub_s.size(); ++i){
			dub_S = somma_triangolare_s(dub_S, dub_s[i]);
		}
		return dub_S;
	}
}

/*
 * gli do una sequenza di tanti l e lui  trova L12; L123; … ;L; restituisce solo con L fisso
 */

inline std::vector <vettore_ushort> somma_triangolare_s (const std::vector <unsigned short> & dub_s, 
														 const unsigned short & dub_S) {
	
	std::vector <vettore_ushort> dub_S_temp (somma_triangolare_s (dub_s));
	std::vector <vettore_ushort> dub_S_tot;
	for (short i = 0 ; i < dub_S_temp.size() ; ++i) {
		if (dub_S_temp[i].back()==dub_S) {
			dub_S_tot.push_back(dub_S_temp[i]);
		}
		else;
	}
	return dub_S_tot;	
}


#endif

