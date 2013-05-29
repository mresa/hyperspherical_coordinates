/*
 accoppiamenti_Spinoriali 
 
 
 
 
 BASE 1 s,sz
 BASE 2 S12 S123 ... S , SZ
 
 
 *  Created by Marco Resa on 13/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef ACC_S_H
#define ACC_S_H

#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "algebra.h"
#include "accoppiamenti_iperangolari.h"
#include "clebsch_s.h"


/*
 * 
 * dato [S], Z, [s] ne faccio lo sviluppo nella base s1,z1,s2,z2 … sA, zA con i relativi coefficienti
 * 
 * 
 */

class Accoppiamenti_Spinoriali {
	
	unsigned short AA;											// A totale
	std::vector < unsigned short > dub_ss;							
	std::vector < unsigned short > dub_SS;							
	short dub_ZZ;												// M totale
	unsigned short degenerazione;								// numero di stati totali
	std::vector < long double > coefficienti;
	std::vector < vettore_short > dub_zz;						// vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	
	
public:
	
	
	Accoppiamenti_Spinoriali (const std::vector < unsigned short >&  dub_s , 
							  const std::vector < unsigned short >& dub_S ,
							  const short&  dub_Z) {
		if (dub_s.size() == 1){
			assert (dub_s==dub_S); //blocca lo script, da corregere?
			AA = dub_s.size();
			dub_ss = dub_s;
			dub_SS = dub_S;
			dub_ZZ = dub_Z;
			degenerazione = 1;
			std::vector <short> dub_zz_temp(1,dub_ZZ);
			dub_zz.push_back(dub_zz_temp);
			coefficienti.push_back(1);
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			if(controlla_s_z(dub_S.back(), dub_Z) && dis_triang_s(dub_S, dub_s)){
				AA = dub_s.size();
				dub_ss = dub_s;
				dub_SS = dub_S;
				dub_ZZ = dub_Z;
				degenerazione = 0;
			
				std::vector <vettore_short> dub_ZZ_temp (costruisci_Z (dub_SS , dub_ZZ));
				for(short i = -dub_ss[0] ; i <= dub_ss[0] ; i+=2){
					for(short j = 0 ; j < dub_ZZ_temp.size() ; ++j){
						std::vector <short> dub_zz_temp (m_combinazioni_di_M(dub_ZZ_temp[j] , i));
						if (controlla_s_z(dub_ss,dub_zz_temp)){
							long double coeff = clebsch_gordan_s (dub_ss, dub_zz_temp, dub_SS, dub_ZZ);
							if (coeff != 0) {
								dub_zz.push_back(dub_zz_temp);
								coefficienti.push_back(coeff);
							}
							else;
						}
						else;
					}
				}
				degenerazione = coefficienti.size();
			}
			else{
				AA = dub_s.size();
				dub_SS = dub_S;
				dub_ZZ = dub_Z;
				dub_ss = dub_s;
				degenerazione = 0;
			}
		}
	}
	
	//distruttore
	~Accoppiamenti_Spinoriali () {
	}
	
	// N
	unsigned short A() const {return AA;}
	
	unsigned short D() const {return degenerazione;}	
	// S totale
	std::vector< unsigned short > D_S() const {return dub_SS;}
	unsigned short D_S(const unsigned short& i) const {
		assert( i < AA - 1  );
		return dub_SS[i];
	}	

	unsigned short D_S_tot() const {return dub_SS.back();}
	
	// M totale
	short D_Z() const {return dub_ZZ;}
	
	
	// accesso agli s
	std::vector< unsigned short > D_s() const {return dub_ss;}
	unsigned short D_s(const unsigned short& i) const {
		assert( i < AA );
		return dub_ss[i];
	}	
		
	// accesso agli m
	std::vector< vettore_short > D_z() const {return dub_zz;}
	std::vector< short > D_z(const unsigned short& i) const {
		assert( i < degenerazione );
		return dub_zz[i];
	}
	short D_z(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < AA );
		return dub_zz[i][j];
	}	
	
	// accesso ai coefficienti
	std::vector< long double > cm() const {
		return coefficienti;
	}
	long double cm(const unsigned short& i) const{
		assert( i < degenerazione );
		return coefficienti[i];
	}
};

/*
 * 
 * dato [S], Z, [s] ne faccio lo sviluppo nella base s1,z1,s2,z2 … sA, zA con i relativi coefficienti, clusterizzando
 * 
 * 
 */
class Accoppiamenti_Spinoriali_cluster {
	
	unsigned short AA;											// A totale
	std::vector < vettore_ushort > dub_ss;							
	std::vector < vettore_ushort > dub_SS_r;							
	std::vector < unsigned short > dub_SS;							
	short dub_ZZ;												// M totale
	unsigned short degenerazione;								// numero di stati totali
	std::vector < long double > coefficienti;
	std::vector < v_vettore_short > dub_zz;						// vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	
	
public:
	
	
	Accoppiamenti_Spinoriali_cluster (const std::vector < vettore_ushort > & dub_s , 
									  const std::vector < vettore_ushort > & dub_S_r,
									  const std::vector < unsigned short > & dub_S ,
									  const short &  dub_Z) {
		if(controlla_s_z(dub_S.back(), dub_Z) && dis_triang_s(dub_S_r, dub_s)  && dis_triang_s(dub_S, dub_S_r)){
			AA = 0;
			for(short spazi = 0 ; spazi < dub_s.size() ; ++spazi){
				AA += dub_s[spazi].size();
			}
			
			dub_ss = dub_s;
			dub_SS_r = dub_S_r;
			dub_SS = dub_S;
			dub_ZZ = dub_Z;
			degenerazione = 0;
			
			
			std::vector <v_vettore_short> dub_zz_temp (costruisci_Z(dub_ss));
			
			for(short i = 0 ; i < dub_zz_temp.size() ; ++i){
				long double coeff = clebsch_gordan_s_cluster (dub_ss, dub_zz_temp[i], dub_SS_r, dub_SS, dub_ZZ);
				if (coeff != 0) {
					dub_zz.push_back(dub_zz_temp[i]);
					coefficienti.push_back(coeff);
				}
				else;
			}
			
			degenerazione = coefficienti.size();		
		}
		else{
			AA = dub_s.size();
			dub_SS = dub_S;
			dub_ZZ = dub_Z;
			dub_ss = dub_s;
			degenerazione = 0;
		}
	}
	
	
	//distruttore
	~Accoppiamenti_Spinoriali_cluster () {
	}
	
	// N
	unsigned short A() const {return AA;}
	
	unsigned short D() const {return degenerazione;}	
	// S totale
	std::vector< vettore_ushort > D_S_v() const {return dub_SS_r;}
	std::vector< unsigned short > D_S_r() const {
		return declusterizza(dub_SS_r);
	}
	std::vector< unsigned short > D_S_t() const {return dub_SS;}
	std::vector< unsigned short > D_S() const {
		std::vector< unsigned short > D_S_temp(declusterizza(dub_SS_r));
		for (short i = 0; i<dub_SS.size(); ++i) {
			D_S_temp.push_back(dub_SS[i]);
		}
		return D_S_temp;
	}
	
	
	unsigned short D_S_tot() const {return dub_SS.back();}
	
	// M totale
	short D_Z() const {return dub_ZZ;}
	
	
	// accesso agli s
	std::vector< vettore_ushort > D_s_v() const {return dub_ss;}
	std::vector< unsigned short > D_s() const {
		return declusterizza(dub_ss);
	}
	
	// accesso agli m
	std::vector< v_vettore_short > D_z_v() const {return dub_zz;}
	std::vector< vettore_short > D_z() const {
		std::vector < vettore_short > dub_ss_temp;
		for(short casi = 0 ; casi < dub_zz.size() ; ++casi){
			dub_ss_temp.push_back(declusterizza(dub_zz[casi]));
		}
		return dub_ss_temp;
	}
	
	// accesso ai coefficienti
	std::vector< long double > cm() const {
		return coefficienti;
	}
	long double cm(const unsigned short& i) const{
		assert( i < degenerazione );
		return coefficienti[i];
	}
};

/*
 * 
 * dati gli stati [s]_a genero i possibili stati [S]_a-1
 * dati gli stati [S]_a-1 [s]_a genero i possibili stati [S]_a-1  [s]_a Z  -> 2a
 *
 * per ogni stato [S],[s]_a, M genero lo sviluppo in [s][s_z] 
 * 
 * eventualmente a S fisso o S,Z fisso.
 */


class Accoppiamenti_Spinoriali_S {
	
	unsigned short AA;											// A totale
	std::vector < unsigned short > dub_ss;						// [s]	
	std::vector < vettore_ushort > dub_SS;						// [S] possibili _i
	std::vector < short > dub_ZZ;								// possibili M per _i -> M ij
	unsigned short degenerazione;								// numero di stati i
	std::vector < unsigned short > degenerazione_z;								// numero di stati i
	std::vector < v_vettore_short > dub_zz;						// per ogni stato ij ho k vettori di decomposizione vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	std::vector < vettore_ldouble > coefficienti;				// per ogni stato ij k ho un numero
	std::vector < vettore_bool > parita_stato;						// vettore dei coefficienti per lo sviluppo nella base
	
	static std::vector < bool > genera_parita_stato_s ( const std::vector < unsigned short > &S){
		std::vector < bool > parita;
		parita.push_back((S.front()/2 + 1)%2);
		return parita;
	}
	
	static std::vector < bool > genera_parita_stato_s (std::vector < vettore_ushort > & S_relativo){
		std::vector < bool > parita;
		for (short i = 0 ; i < S_relativo.size() ; ++i) {
			parita.push_back((S_relativo[i].front()/2 + 1 )%2 );
		}
		return parita;
	}
	

public:
	
	
	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1){
			dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			dub_SS_temp = (somma_triangolare_s(dub_s));
		}
		
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_ZZ_temp);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}

	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s, 
								const unsigned short & dub_S) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1 ){
			if (dub_s.back() == dub_S) dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			dub_SS_temp = (somma_triangolare_s(dub_s,dub_S)); //somma triangolare a S fisso
		}
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_ZZ_temp);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}
	
	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s, 
								const unsigned short & dub_S, 
								const short & dub_Z) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1 ){
			if (dub_s.back() == dub_S) dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			dub_SS_temp = (somma_triangolare_s(dub_s,dub_S));
		}
		
		
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			if(controlla_s_z(dub_SS_temp[i].back(), dub_Z)){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_Z);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}

	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s, 
								const unsigned short & dub_S, 
								const short & dub_Z,
								const unsigned short parita) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1 ){
			if (dub_s.back() == dub_S) dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			dub_SS_temp = (somma_triangolare_s(dub_s,dub_S));
		}
		
		
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			if(controlla_s_z(dub_SS_temp[i].back(), dub_Z)){	
				if (dub_SS_temp[i].front()/2. == parita) {
					dub_SS.push_back(dub_SS_temp[i]);
					parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
					dub_ZZ.push_back(dub_Z);
				}
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());				
		}
	}
	
	//cluster	
	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s,
								const std::vector < unsigned short > & CLUSTER) {
		
		if (CLUSTER.size()==1) {
			AA = dub_s.size();
			dub_ss = dub_s;
			std::vector < vettore_ushort > dub_SS_temp;
			if (dub_s.size() == 1){
				dub_SS_temp.push_back(dub_ss) ;
			}
			else{
				assert(dub_s.size() > 1); // almeno 2 o non funziona
				dub_SS_temp = (somma_triangolare_s(dub_s));
			}
			
			
			for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
				for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
					dub_SS.push_back(dub_SS_temp[i]);
					parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
					dub_ZZ.push_back(dub_ZZ_temp);
				}
			}
			
			degenerazione = dub_SS.size();			
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
		}
		else{
			AA = dub_ss.size();
			dub_ss=dub_s;
			std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
			unsigned short j_temp = 0;
			for (short i = 0 ; i<CLUSTER.size(); ++i) {
				vettore_ushort dub_s_cl_t;
				for (short j = 0 ; j<CLUSTER[i]; ++j) {
					dub_s_cl_t.push_back(dub_s[j+j_temp]);
				}
				dub_s_cluster.push_back(dub_s_cl_t);
				j_temp += dub_s_cl_t.size();
			}
			
			
			std::vector < v_vettore_ushort > dub_SS_rel_temp;
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
			
			
			for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
				vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
				
				std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
				dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
				dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
			}
			
			for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
				dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
			}
			
			// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
			
			std::vector < vettore_ushort > dub_SS_ultimi;
			for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
				
				vettore_ushort dub_SS_ultimi_t;
				for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
					dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
				}
				dub_SS_ultimi.push_back(dub_SS_ultimi_t);
			}
			
			// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
			
			std::vector < v_vettore_ushort > dub_SS_tot;
			//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
			
			for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
				dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
			}
			
			std::vector < vettore_ushort > dub_SS_tot_t;
			std::vector < v_vettore_ushort > dub_SS_rel_t;
			std::vector < vettore_ushort > dub_SS_ultimi_t;
			
			for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
				for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
					for(short dub_ZZ_temp = - dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp <= dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp += 2){	
						dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
						dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
						dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
						dub_ZZ.push_back(dub_ZZ_temp);
					}
				}
			}
			degenerazione = dub_SS_tot_t.size();			
			for (short i = 0 ; i < degenerazione; ++i) {			
				dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
				for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
					dub_SS[i].push_back(dub_SS_tot_t[i][j]);
				}
			}
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
				parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
			
		}
					
	}
	
	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s,
								const std::vector < unsigned short > & CLUSTER,
								const unsigned short & dub_S) {
		if (CLUSTER.size()==1) {
			AA = dub_s.size();
			dub_ss = dub_s;
			std::vector < vettore_ushort > dub_SS_temp;
			if (dub_s.size() == 1 ){
				if (dub_s.back() == dub_S) dub_SS_temp.push_back(dub_ss) ;
			}
			else{
				assert(dub_s.size() > 1); // almeno 2 o non funziona
				dub_SS_temp = (somma_triangolare_s(dub_s,dub_S)); //somma triangolare a S fisso
			}
			
			for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
				for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
					dub_SS.push_back(dub_SS_temp[i]);
					parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
					dub_ZZ.push_back(dub_ZZ_temp);
				}
			}
			
			degenerazione = dub_SS.size();			
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
		}
		else {
			AA=dub_s.size();
			dub_ss=dub_s;
			std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
			unsigned short j_temp = 0;
			for (short i = 0 ; i<CLUSTER.size(); ++i) {
				vettore_ushort dub_s_cl_t;
				for (short j = 0 ; j<CLUSTER[i]; ++j) {
					dub_s_cl_t.push_back(dub_s[j+j_temp]);
				}
				dub_s_cluster.push_back(dub_s_cl_t);
				j_temp += dub_s_cl_t.size();
			}
			
			std::vector < v_vettore_ushort > dub_SS_rel_temp;
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
			
			
			for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
				vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
				
				std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
				dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
				dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
			}
			
			for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
				dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
			}
			
			// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
			
			std::vector < vettore_ushort > dub_SS_ultimi;
			for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
				
				vettore_ushort dub_SS_ultimi_t;
				for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
					dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
				}
				dub_SS_ultimi.push_back(dub_SS_ultimi_t);
			}
			
			// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
			
			std::vector < v_vettore_ushort > dub_SS_tot;
			//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
			
			for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
				dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
			}
			
			std::vector < vettore_ushort > dub_SS_tot_t;
			std::vector < v_vettore_ushort > dub_SS_rel_t;
			std::vector < vettore_ushort > dub_SS_ultimi_t;
			
			for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
				for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
					if(dub_SS_tot[caso][caso1].back() == dub_S){
						for(short dub_ZZ_temp = - dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp <= dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp += 2){	
							dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
							dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
							dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
							dub_ZZ.push_back(dub_ZZ_temp);
						}
					}
				}
			}
			degenerazione = dub_SS_tot_t.size();			
			for (short i = 0 ; i < degenerazione; ++i) {			
				dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
				for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
					dub_SS[i].push_back(dub_SS_tot_t[i][j]);
				}
			}
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
				parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
		}
	}
	
	Accoppiamenti_Spinoriali_S (const std::vector < unsigned short >&  dub_s,
								const std::vector < unsigned short > & CLUSTER, 
								const unsigned short & dub_S, 
								const short & dub_Z) {
		if (CLUSTER.size()==1) {
			AA = dub_s.size();
			dub_ss = dub_s;
			std::vector < vettore_ushort > dub_SS_temp;
			if (dub_s.size() == 1 ){
				if (dub_s.back() == dub_S) dub_SS_temp.push_back(dub_ss) ;
			}
			else{
				assert(dub_s.size() > 1); // almeno 2 o non funziona
				dub_SS_temp = (somma_triangolare_s(dub_s,dub_S));
			}
			
			
			
			for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
				if(controlla_s_z(dub_SS_temp[i].back(), dub_Z)){	
					dub_SS.push_back(dub_SS_temp[i]);
					parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
					dub_ZZ.push_back(dub_Z);
				}
			}
			
			degenerazione = dub_SS.size();			
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
		}
		else{
			AA=dub_s.size();
			dub_ss=dub_s;
			std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
			unsigned short j_temp = 0;
			for (short i = 0 ; i<CLUSTER.size(); ++i) {
				vettore_ushort dub_s_cl_t;
				for (short j = 0 ; j<CLUSTER[i]; ++j) {
					dub_s_cl_t.push_back(dub_s[j+j_temp]);
				}
				dub_s_cluster.push_back(dub_s_cl_t);
				j_temp += dub_s_cl_t.size();
			}
			
			std::vector < v_vettore_ushort > dub_SS_rel_temp;
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
			
			
			for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
				vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
				
				std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
				dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
				dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
			}
			
			for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
				dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
			}
			
			// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
			
			std::vector < vettore_ushort > dub_SS_ultimi;
			for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
				
				vettore_ushort dub_SS_ultimi_t;
				for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
					dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
				}
				dub_SS_ultimi.push_back(dub_SS_ultimi_t);
			}
			
			// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
			
			std::vector < v_vettore_ushort > dub_SS_tot;
			//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
			
			for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
				dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
			}
			
			std::vector < vettore_ushort > dub_SS_tot_t;
			std::vector < v_vettore_ushort > dub_SS_rel_t;
			std::vector < vettore_ushort > dub_SS_ultimi_t;
			
			for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
				for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
					if(dub_SS_tot[caso][caso1].back() == dub_S){
						dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
						dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
						dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
						dub_ZZ.push_back(dub_Z);
					}
				}
			}
			degenerazione = dub_SS_tot_t.size();			
			for (short i = 0 ; i < degenerazione; ++i) {			
				dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
				for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
					dub_SS[i].push_back(dub_SS_tot_t[i][j]);
				}
			}
			
			for(short i = 0 ; i < degenerazione ; ++i){
				Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
				parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
				dub_zz.push_back(ACC.D_z());
				coefficienti.push_back(ACC.cm());
				degenerazione_z.push_back(ACC.D());
			}
		}
	}
	
	//distruttore
	~Accoppiamenti_Spinoriali_S () {
	}
	
	
	inline unsigned short A() const {return AA;}
	
	inline unsigned short D() const {return degenerazione;}	
	
	inline std::vector <unsigned short> Dz() const {return degenerazione_z;}	
	inline unsigned short Dz(const unsigned short& i) const {
		assert( i < degenerazione );
		return degenerazione_z[i];
	}	

	// S totale
	inline std::vector< vettore_ushort > D_S() const {return dub_SS;}
	inline std::vector< unsigned short > D_S(const unsigned short & i) const {
		assert( i < degenerazione  );
		return dub_SS[i];
	}
	
	std::vector < vettore_bool > P_S() const {return parita_stato;}
	std::vector< bool > P_S(const unsigned short& i) const {
		assert( i < degenerazione );
		return parita_stato[i];
	}
	
	inline unsigned short D_S(const unsigned short& i, const unsigned short& j) const {
		assert( i < degenerazione  );
		assert( j < AA - 1  );
		return dub_SS[i][j];
	}	

	inline std::vector< short > D_Z() const {return dub_ZZ;}
	inline short D_Z(const unsigned short & i) const {
		assert( i < degenerazione );
		return dub_ZZ[i];
	}
	
	// accesso agli s
	inline std::vector< unsigned short > D_s() const {return dub_ss;}
	inline unsigned short D_s(const unsigned short& i) const {
		assert( i < AA );
		return dub_ss[i];
	}	
	

	// accesso agli m
	inline std::vector< v_vettore_short > D_z() const {return dub_zz;}
	inline std::vector< vettore_short > D_z(const unsigned short& i) const {
		assert( i < degenerazione );
		return dub_zz[i];
	}
	inline std::vector< short > D_z(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < dub_zz[i].size()  );
		return dub_zz[i][j];
	}
	inline short D_z(const unsigned short& i,const unsigned short& j, const unsigned short& k) const {
		assert( i < degenerazione );
		assert( j < dub_zz[i].size() );
		assert( k < AA );
		return dub_zz[i][j][k];
	}	

	// accesso ai coefficienti
	inline std::vector< vettore_ldouble > cm() const {return coefficienti;}
	inline std::vector< long double > cm(const unsigned short& i) const {
		assert( i < degenerazione );
		return coefficienti[i];
	}
	inline long double cm(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < coefficienti[i].size());
		return coefficienti[i][j];
	}
};


/*
 * 
 * Come prima, sviluppo fino a S_max
 *
 * Eventualmente a S_z Fisso (utile per isospin)
 */


class Accoppiamenti_Spinoriali_S_max {
	
	unsigned short AA;											// A totale
	std::vector < unsigned short > dub_ss;						// [s]	
	std::vector < vettore_ushort > dub_SS;						// [S] possibili _i
	std::vector < short > dub_ZZ;								// possibili M per _i -> M ij
	unsigned short degenerazione;								// numero di stati i
	std::vector < unsigned short > degenerazione_z;				// numero di stati i
	std::vector < v_vettore_short > dub_zz;						// per ogni stato ij ho k vettori di decomposizione vettore dei valori [m] : tra i (2*l1+1)(2*l2+1) vettori a 2 componenti possibili scelgo quelli tali che M = m1 + m2
	std::vector < vettore_ldouble > coefficienti;				// per ogni stato ij k ho un numero
	std::vector < vettore_bool > parita_stato;						// vettore dei coefficienti per lo sviluppo nella base

	
	static std::vector < bool > genera_parita_stato_s ( const std::vector < unsigned short > &S){
		std::vector < bool > parita;
		parita.push_back((S.front()/2 + 1)%2);
		return parita;
	}
	
	static std::vector < bool > genera_parita_stato_s (std::vector < vettore_ushort > & S_relativo){
		std::vector < bool > parita;
		for (short i = 0 ; i < S_relativo.size() ; ++i) {
			parita.push_back((S_relativo[i].front()/2 + 1 )%2 );
		}
		return parita;
	}
	
public:
	
	
/*	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >& dub_s) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1){
			dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			dub_SS_temp = (somma_triangolare_s(dub_s));
		}
		
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_ZZ_temp);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >& dub_s, 
									const unsigned short & dub_S_max) {
		
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1 ){
			if (dub_s.back() <= dub_S_max) dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			unsigned short dub_S_parita=AA%2; 
			//std::vector < vettore_ushort > dub_SS_temp;// = (somma_triangolare_s(dub_s,dub_S_parita)); //somma triangolare a S fisso

			dub_SS_temp = (somma_triangolare_s(dub_s,dub_S_parita)); //somma triangolare a S fisso

						
			for (short i = dub_S_parita + 2; i<= dub_S_max; i+=2) {
				std::vector < vettore_ushort > dub_SS_temp_t = (somma_triangolare_s(dub_s,i));
				for (short j = 0; j < dub_SS_temp_t.size(); j++) {
					dub_SS_temp.push_back(dub_SS_temp_t[j]);
				}
			}
		}
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			for(short dub_ZZ_temp = - dub_SS_temp[i].back() ; dub_ZZ_temp <= dub_SS_temp[i].back() ; dub_ZZ_temp += 2){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_ZZ_temp);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s, 
									const unsigned short & dub_S_max, 
									const short & dub_Z) {
		AA = dub_s.size();
		dub_ss = dub_s;
		std::vector < vettore_ushort > dub_SS_temp;
		if (dub_s.size() == 1 ){
			if (dub_s.back() <= dub_S_max) dub_SS_temp.push_back(dub_ss) ;
		}
		else{
			assert(dub_s.size() > 1); // almeno 2 o non funziona
			unsigned short dub_S_parita=AA%2; 
			dub_SS_temp = (somma_triangolare_s(dub_s,dub_S_parita)); //somma triangolare a S fisso
			for (short i = dub_S_parita + 2; i<= dub_S_max; i+=2) {
				std::vector < vettore_ushort > dub_SS_temp_t = (somma_triangolare_s(dub_s,i));
				for (short j = 0; j < dub_SS_temp_t.size(); j++) {
					dub_SS_temp.push_back(dub_SS_temp_t[j]);
				}
			}
		}
		
		
		
		for(short i = 0 ; i < dub_SS_temp.size() ; ++i){
			if(controlla_s_z(dub_SS_temp[i].back(), dub_Z)){	
				dub_SS.push_back(dub_SS_temp[i]);
				parita_stato.push_back(genera_parita_stato_s(dub_SS.back()));
				dub_ZZ.push_back(dub_Z);
			}
		}
		
		degenerazione = dub_SS.size();			
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali ACC(dub_ss,dub_SS[i],dub_ZZ[i]);
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
	}

	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER) {
		
		dub_ss=dub_s;
		std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
		unsigned short j_temp = 0;
		for (short i = 0 ; i<CLUSTER.size(); ++i) {
			vettore_ushort dub_s_cl_t;
			for (short j = 0 ; j<CLUSTER[i]; ++j) {
				dub_s_cl_t.push_back(dub_s[j+j_temp]);
			}
			dub_s_cluster.push_back(dub_s_cl_t);
			j_temp += dub_s_cl_t.size();
		}
		
		std::vector < v_vettore_ushort > dub_SS_rel_temp;
		std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
		
		
		for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
			vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
			
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
			dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
			dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
		}
		
		for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
			dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
		}
		
		// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
		
		std::vector < vettore_ushort > dub_SS_ultimi;
		for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
			
			vettore_ushort dub_SS_ultimi_t;
			for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
				dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
			}
			dub_SS_ultimi.push_back(dub_SS_ultimi_t);
		}
		
		// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
		
		std::vector < v_vettore_ushort > dub_SS_tot;
		//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
		
		for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
			dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
		}
		
		std::vector < vettore_ushort > dub_SS_tot_t;
		std::vector < v_vettore_ushort > dub_SS_rel_t;
		std::vector < vettore_ushort > dub_SS_ultimi_t;
		
		for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
			for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
				for(short dub_ZZ_temp = - dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp <= dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp += 2){	
					dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
					dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
					dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
					dub_ZZ.push_back(dub_ZZ_temp);
				}
			}
		}
		degenerazione = dub_SS_tot_t.size();			
		for (short i = 0 ; i < degenerazione; ++i) {			
			dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
			for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
				dub_SS[i].push_back(dub_SS_tot_t[i][j]);
			}
		}
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
			parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
		
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER,
									const unsigned short & dub_S_max) {
		
		dub_ss=dub_s;
		std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
		unsigned short j_temp = 0;
		for (short i = 0 ; i<CLUSTER.size(); ++i) {
			vettore_ushort dub_s_cl_t;
			for (short j = 0 ; j<CLUSTER[i]; ++j) {
				dub_s_cl_t.push_back(dub_s[j+j_temp]);
			}
			dub_s_cluster.push_back(dub_s_cl_t);
			j_temp += dub_s_cl_t.size();
		}
		
		std::vector < v_vettore_ushort > dub_SS_rel_temp;
		std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
		
		
		for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
			vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
			
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
			dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
			dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
		}
		
		for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
			dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
		}
		
		// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
		
		std::vector < vettore_ushort > dub_SS_ultimi;
		for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
			
			vettore_ushort dub_SS_ultimi_t;
			for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
				dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
			}
			dub_SS_ultimi.push_back(dub_SS_ultimi_t);
		}
		
		// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
		
		std::vector < v_vettore_ushort > dub_SS_tot;
		//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
		
		for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
			dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
		}
		
		std::vector < vettore_ushort > dub_SS_tot_t;
		std::vector < v_vettore_ushort > dub_SS_rel_t;
		std::vector < vettore_ushort > dub_SS_ultimi_t;
		
		for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
			for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
				if(dub_SS_tot[caso][caso1].back() <= dub_S_max){
					for(short dub_ZZ_temp = - dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp <= dub_SS_tot[caso][caso1].back() ; dub_ZZ_temp += 2){	
						dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
						dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
						dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
						dub_ZZ.push_back(dub_ZZ_temp);
					}
				}
			}
		}
		degenerazione = dub_SS_tot_t.size();			
		for (short i = 0 ; i < degenerazione; ++i) {			
			dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
			for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
				dub_SS[i].push_back(dub_SS_tot_t[i][j]);
			}
		}
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
			parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
		
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER, 
									const unsigned short & dub_S_max, 
									const short & dub_Z) {
		
		dub_ss=dub_s;
		std::vector < vettore_ushort > dub_s_cluster; //vettore degli s delle partizioni
		unsigned short j_temp = 0;
		for (short i = 0 ; i<CLUSTER.size(); ++i) {
			vettore_ushort dub_s_cl_t;
			for (short j = 0 ; j<CLUSTER[i]; ++j) {
				dub_s_cl_t.push_back(dub_s[j+j_temp]);
			}
			dub_s_cluster.push_back(dub_s_cl_t);
			j_temp += dub_s_cl_t.size();
		}
		
		std::vector < v_vettore_ushort > dub_SS_rel_temp;
		std::vector < vettore_ushort > dub_SS_rel_temp_inizio(somma_triangolare_s(dub_s_cluster[0]));
		
		
		for (short casi = 0 ; casi< dub_SS_rel_temp_inizio.size(); ++casi) {
			vettore_ushort dub_SS_rel_temp_inizio_temp(dub_SS_rel_temp_inizio[casi]);
			
			std::vector < vettore_ushort > dub_SS_rel_temp_inizio_temp_t;
			dub_SS_rel_temp_inizio_temp_t.push_back(dub_SS_rel_temp_inizio_temp);
			dub_SS_rel_temp.push_back(dub_SS_rel_temp_inizio_temp_t);
		}
		
		for (short spazio = 1 ; spazio<CLUSTER.size(); ++spazio) {
			dub_SS_rel_temp = somma_diretta_M(dub_SS_rel_temp,somma_triangolare_s(dub_s_cluster[spazio]));
		}
		
		// ho popolato dub_SS_rel_temp un  vettore sui casi, uno sugli spazi, uno sui numeri
		
		std::vector < vettore_ushort > dub_SS_ultimi;
		for(short caso = 0 ; caso < dub_SS_rel_temp.size() ; ++caso){
			
			vettore_ushort dub_SS_ultimi_t;
			for(short spazio = 0 ; spazio < dub_SS_rel_temp[caso].size() ; ++spazio){
				dub_SS_ultimi_t.push_back(dub_SS_rel_temp[caso][spazio].back());
			}
			dub_SS_ultimi.push_back(dub_SS_ultimi_t);
		}
		
		// ho popolato gli ultimi cerco l'accoppiamento tra i cluster
		
		std::vector < v_vettore_ushort > dub_SS_tot;
		//1 sui casi di dub_SS_ultimi, 2 sui casi di ora, 3 sui numeri
		
		for (short caso = 0; caso < dub_SS_ultimi.size(); ++caso) {
			dub_SS_tot.push_back(somma_triangolare_s(dub_SS_ultimi[caso]));
		}
		
		std::vector < vettore_ushort > dub_SS_tot_t;
		std::vector < v_vettore_ushort > dub_SS_rel_t;
		std::vector < vettore_ushort > dub_SS_ultimi_t;
		
		for(short caso = 0 ; caso < dub_SS_ultimi.size() ; ++caso){
			for (short caso1 = 0; caso1 < dub_SS_tot[caso].size(); ++caso1) {
				if(dub_SS_tot[caso][caso1].back() <= dub_S_max && controlla_s_z(dub_SS_tot[caso][caso1].back(), dub_Z)){
					dub_SS_rel_t.push_back(dub_SS_rel_temp[caso]);
					dub_SS_ultimi_t.push_back(dub_SS_ultimi[caso]);
					dub_SS_tot_t.push_back(dub_SS_tot[caso][caso1]);
					dub_ZZ.push_back(dub_Z);
				}
			}
		}
		degenerazione = dub_SS_tot_t.size();			
		for (short i = 0 ; i < degenerazione; ++i) {			
			dub_SS.push_back(declusterizza(dub_SS_rel_t[i]));
			for (short j = 0 ; j < dub_SS_tot_t[i].size(); ++j) {
				dub_SS[i].push_back(dub_SS_tot_t[i][j]);
			}
		}
		
		for(short i = 0 ; i < degenerazione ; ++i){
			Accoppiamenti_Spinoriali_cluster ACC(dub_s_cluster,dub_SS_rel_t[i],dub_SS_tot_t[i],dub_ZZ[i]);
			parita_stato.push_back(genera_parita_stato_s( dub_SS_rel_t[i]));
			dub_zz.push_back(ACC.D_z());
			coefficienti.push_back(ACC.cm());
			degenerazione_z.push_back(ACC.D());
		}
		
	}
*/	
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >& dub_s) {
		AA = dub_s.size();
		dub_ss = dub_s;
		Accoppiamenti_Spinoriali_S ACC(dub_ss);
		dub_SS = ACC.D_S();
		dub_ZZ = ACC.D_Z();
		degenerazione = ACC.D();
		degenerazione_z = ACC.Dz();
		dub_zz = ACC.D_z();
		coefficienti = ACC.cm();
		parita_stato = ACC.P_S();
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >& dub_s, 
									const unsigned short & dub_S_max) {
		AA = dub_s.size();
		dub_ss = dub_s;
		bool dub_S_min = AA%2;
		degenerazione = 0;
		for (unsigned short dub_S_temp = dub_S_min; dub_S_temp <= dub_S_max; dub_S_temp+=2) {
			Accoppiamenti_Spinoriali_S ACC(dub_ss,dub_S_temp);
			degenerazione += ACC.D();
			for (unsigned long i = 0; i<ACC.D(); ++i) {
				dub_SS.push_back(ACC.D_S(i));
				dub_ZZ.push_back(ACC.D_Z(i));
				degenerazione_z.push_back(ACC.Dz(i));
				dub_zz.push_back(ACC.D_z(i));
				coefficienti.push_back(ACC.cm(i));
				parita_stato.push_back(ACC.P_S(i));			
			}
		}
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s, 
									const unsigned short & dub_S_max, 
									const short & dub_Z) {
		AA = dub_s.size();
		dub_ss = dub_s;
		degenerazione = 0;
		bool dub_S_min = AA%2;
		for (unsigned short dub_S_temp = dub_S_min; dub_S_temp <= dub_S_max; dub_S_temp+=2) {
			if (controlla_s_z(dub_S_temp, dub_Z)) {
				Accoppiamenti_Spinoriali_S ACC(dub_ss,dub_S_temp,dub_Z);
				degenerazione += ACC.D();
				for (unsigned long i = 0; i<ACC.D(); ++i) {
					dub_SS.push_back(ACC.D_S(i));
					dub_ZZ.push_back(ACC.D_Z(i));
					degenerazione_z.push_back(ACC.Dz(i));
					dub_zz.push_back(ACC.D_z(i));
					coefficienti.push_back(ACC.cm(i));
					parita_stato.push_back(ACC.P_S(i));			
				}				
			}
		}
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER) {
		AA = dub_s.size();
		dub_ss = dub_s;
		Accoppiamenti_Spinoriali_S ACC(dub_ss,CLUSTER);
		dub_SS = ACC.D_S();
		dub_ZZ = ACC.D_Z();
		degenerazione = ACC.D();
		degenerazione_z = ACC.Dz();
		dub_zz = ACC.D_z();
		coefficienti = ACC.cm();
		parita_stato = ACC.P_S();
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER,
									const unsigned short & dub_S_max) {
		AA = dub_s.size();
		dub_ss = dub_s;
		degenerazione = 0;
		bool dub_S_min = AA%2;
		for (unsigned short dub_S_temp = dub_S_min; dub_S_temp <= dub_S_max; dub_S_temp+=2) {
			Accoppiamenti_Spinoriali_S ACC(dub_ss,CLUSTER,dub_S_temp);
			degenerazione += ACC.D();
			for (unsigned long i = 0; i<ACC.D(); ++i) {
				dub_SS.push_back(ACC.D_S(i));
				dub_ZZ.push_back(ACC.D_Z(i));
				degenerazione_z.push_back(ACC.Dz(i));
				dub_zz.push_back(ACC.D_z(i));
				coefficienti.push_back(ACC.cm(i));
				parita_stato.push_back(ACC.P_S(i));			
			}
		}
	}
	
	Accoppiamenti_Spinoriali_S_max (const std::vector < unsigned short >&  dub_s,
									const std::vector < unsigned short > & CLUSTER, 
									const unsigned short & dub_S_max, 
									const short & dub_Z) {
		AA = dub_s.size();
		dub_ss = dub_s;
		degenerazione = 0;
		bool dub_S_min = AA%2;
		for (unsigned short dub_S_temp = dub_S_min; dub_S_temp <= dub_S_max; dub_S_temp+=2) {
			if (controlla_s_z(dub_S_temp, dub_Z)) {
				Accoppiamenti_Spinoriali_S ACC(dub_ss,CLUSTER,dub_S_temp,dub_Z);
				degenerazione += ACC.D();
				for (unsigned long i = 0; i<ACC.D(); ++i) {
					dub_SS.push_back(ACC.D_S(i));
					dub_ZZ.push_back(ACC.D_Z(i));
					degenerazione_z.push_back(ACC.Dz(i));
					dub_zz.push_back(ACC.D_z(i));
					coefficienti.push_back(ACC.cm(i));
					parita_stato.push_back(ACC.P_S(i));			
				}				
			}
		}
	}
	
	//distruttore
	~Accoppiamenti_Spinoriali_S_max () {
	}
	
	inline unsigned short A() const {return AA;}
	
	// S totale
	inline std::vector< vettore_ushort > D_S() const {return dub_SS;}
	inline std::vector< unsigned short > D_S(const unsigned short & i) const {
		assert( i < degenerazione  );
		return dub_SS[i];
	}
	inline unsigned short D_S(const unsigned short& i, const unsigned short& j) const {
		assert( i < degenerazione  );
		assert( j < AA - 1  );
		return dub_SS[i][j];
	}	
	std::vector < vettore_bool > P_S() const {return parita_stato;}
	std::vector< bool > P_S(const unsigned short& i) const {
		assert( i < degenerazione );
		return parita_stato[i];
	}
	
	inline std::vector< short > D_Z() const {return dub_ZZ;}
	inline short D_Z(const unsigned short & i) const {
		assert( i < degenerazione );
		return dub_ZZ[i];
	}
	
	// accesso agli s
	inline std::vector< unsigned short > D_s() const {return dub_ss;}
	inline unsigned short D_s(const unsigned short& i) const {
		assert( i < AA );
		return dub_ss[i];
	}	
	
	inline unsigned short D() const {return degenerazione;}	
	inline std::vector <unsigned short> Dz() const {return degenerazione_z;}	
	inline unsigned short Dz(const unsigned short& i) const {
		assert( i < degenerazione );
		return degenerazione_z[i];
	}	
	
	// accesso agli m
	inline std::vector< v_vettore_short > D_z() const {return dub_zz;}
	inline std::vector< vettore_short > D_z(const unsigned short& i) const {
		assert( i < degenerazione );
		return dub_zz[i];
	}
	inline std::vector< short > D_z(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < dub_zz[i].size()  );
		return dub_zz[i][j];
	}
	inline short D_z(const unsigned short& i,const unsigned short& j, const unsigned short& k) const {
		assert( i < degenerazione );
		assert( j < dub_zz[i].size() );
		assert( k < AA );
		return dub_zz[i][j][k];
	}	
	// accesso ai coefficienti
	inline std::vector< vettore_ldouble > cm() const {return coefficienti;}
	inline std::vector< long double > cm(const unsigned short& i) const {
		assert( i < degenerazione );
		return coefficienti[i];
	}
	inline long double cm(const unsigned short& i,const unsigned short& j) const {
		assert( i < degenerazione );
		assert( j < coefficienti[i].size());
		return coefficienti[i][j];
	}
};

#endif

