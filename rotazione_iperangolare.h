/*
 *  rotazione_iperangolare
 *  Tesi
 *
 *  Created by Marco Resa on 29/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ROT_IPER_H
#define ROT_IPER_H

#include "accoppiamenti_iperangolari.h"
#include "tipi_definiti.h"
#include "vettore_rotazione_cinematica.h"
#include "matrice_trasformazione_jacobi.h"
#include "hyper.h"
#include "hyperintegrazione.h"
#include "integrale_angolare_n.h"
#include <omp.h> //calcolo parallelo OMP
#include <boost/timer.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <fstream>
#include "io_sparse.h"
#include "lista_integrazioni.h"



inline std::vector<long double> C_ij(const Accoppiamenti_Iperangolari & TEMP,
									 const unsigned short & indice,
									 const std::vector< vettore_ldouble > & ROT_CIN) {
	unsigned short coppie = ROT_CIN.size();	//sono il numero di coppie, non A
	std::vector<long double> C_ij_t;
	
	if(TEMP.D() > 0){
		hyper H_k_0 (TEMP.n(indice),TEMP.l(indice),TEMP.m(indice,0));
		for(short i = 0; i < coppie ; ++i){
			C_ij_t.push_back(H_k_0.RIPELLE(ROT_CIN[i]));
		} 
	}
	return C_ij_t;
}



MATRICI_SPARSE G_M(const Accoppiamenti_Iperangolari & TEMP, // k'
				   const unsigned short & indice, // k'
				   const HyperInt & G,
				   Integrazioni_3_P * LISTA,
				   const Accoppiamenti_Iperangolari_K_max & K,
				   const unsigned short & dimensione_calcolata) {
	
	assert(indice < TEMP.D());
	
	MATRICI_SPARSE M(K.D(),K.D());
	if(TEMP.D()>0){
		hyper H_k_0 (TEMP.n(indice),TEMP.l(indice),TEMP.m(indice,0));
		for(short i = 0 ; i < K.D() ; ++i){
			hyper H1_0 (K.n(i),K.l(i),K.m(i,0));
			for(short j = i ; j < K.D() ; ++j){				
				if( i >= dimensione_calcolata || j >= dimensione_calcolata){
					hyper H2_0 (K.n(j),K.l(j),K.m(j,0));
					if (dis_triang_H_l(H2_0,H_k_0,H1_0) && dis_triang (K.L(j), TEMP.L(indice), K.L(i))) {
						//interroga la lista ed ottiene il risultato
						//long double hyper_int_temp = hyper_int_ridotto_l(H2_0, H_k_0, H1_0);
						long double hyper_int_temp = 1.;
						
						for (short spindex = 0 ; (spindex < TEMP.N() && hyper_int_temp != 0); ++spindex) {
							hyper_int_temp *= LISTA->Richiesta_L(spindex,H2_0, H_k_0, H1_0) ;
						}
						
						if(hyper_int_temp != 0){
							for (short spindex = 1 ; (spindex < TEMP.N() && hyper_int_temp != 0); ++spindex) {
								hyper_int_temp *= LISTA->Richiesta_P(G,spindex,H2_0,H_k_0,H1_0);
							}	
							if(hyper_int_temp != 0){
								long double temp_m;
								
								bool almk_calcolato;
								temp_m = LISTA->Richiesta_A_LMK(almk_calcolato,H2_0.l(),H_k_0.l(),H1_0.l());
								
								if (almk_calcolato) {
									hyper_int_temp *= temp_m;
									if(std::abs(hyper_int_temp) > errore_matrici){
										M(i,j) = hyper_int_temp;
										M(j,i) = hyper_int_temp;
									}
								}
								else {
#pragma omp parallel for reduction(+:temp_m) schedule( guided )
									for(short d_k = 0; d_k < TEMP.DM(indice); ++d_k){
										hyper H_k (TEMP.n(indice),TEMP.l(indice),TEMP.m(indice,d_k));
										
										long double C1 = TEMP.clmk(indice,d_k);
										
										long double temp1 = 0.;
										if (C1 != 0.){
#pragma omp parallel for reduction(+:temp1) schedule( guided )
											for(short d1 = 0; d1 < K.DM(i); ++d1){
												long double temp2 = 0.;
												hyper H1 (K.n(i),K.l(i),K.m(i,d1));
#pragma omp parallel for reduction(+:temp2) schedule( guided )
												for(short d2 = 0; d2 < K.DM(j); ++d2){
													hyper H2 (K.n(j),K.l(j),K.m(j,d2));
													long double val = 0.;
													if (dis_triang_H_m(H2,H_k,H1)){
														long double integrale = 1.;
#pragma omp critical
														for (short spindex = 0; spindex< TEMP.N() ; ++spindex) {
															integrale *= LISTA->Richiesta_M(spindex,H2, H_k, H1);
															if (integrale == 0) spindex = TEMP.N();
														}
														val = integrale;
														
													}
													temp2 += K.cm(j,d2) * val;
												}
												temp1 += K.cm(i,d1) * temp2;
											}
											temp1 *= C1;
										}
										temp_m += temp1;
									}
									
									LISTA->Richiesta_A_LMK(H2_0.l(),H_k_0.l(),H1_0.l(),temp_m);
									hyper_int_temp *= temp_m;
									if(std::abs(hyper_int_temp) > errore_matrici){
										M(i,j) = hyper_int_temp;
										M(j,i) = hyper_int_temp;
									}
								}
							}
						}
					}
				}
			}
		}
	}
	return M;
}

MATRICI_SPARSE G_M(const Accoppiamenti_Iperangolari & TEMP, // k'
				   const unsigned short & indice, // k'
				   const HyperInt & G,
				   Integrazioni_3_P * LISTA,
				   const Accoppiamenti_Iperangolari_K_max & K) {
	return G_M(TEMP,indice,G,LISTA,K,0);
}

inline V_MATRICI_SPARSE G_M_ij(const unsigned short & n,
							   const Accoppiamenti_Iperangolari & TEMP,
							   const HyperInt & G,
							   Integrazioni_3_P * LISTA,
							   const Accoppiamenti_Iperangolari_K_max & K,
							   const std::vector< vettore_ldouble > & ROT_CIN){
	unsigned short N = TEMP.N();
	V_MATRICI_SPARSE SUPERMATRICE;
	hyper H0 (N,n); // PB
	long double C = (H0.RIPELLE( N - 1, 0. ));
	
	MATRICI_SPARSE G_M_bis_temp_0;
	
	std::fstream f;
	f.precision(std::numeric_limits<long double>::digits10);
	std::string s = "./matrici_G/G_M_base_" + K.fileD() + "_Kg_" + scrivi(TEMP.K())+".txt";
	bool costruisci_file = false;
	f.open(s.c_str(),std::ios::in);
	if (f) {
		std::cerr << "Leggo " << s << "  :";
		unsigned long k=0;
		while (f >> boost::numeric::ublas::io::sparse(G_M_bis_temp_0)){
			std::cerr << " " << k << " ";
			if (k==0) {
				std::vector<long double> C_ij_bis_temp_0 (C_ij(TEMP,0,ROT_CIN)); //uno per coppia
				for(short i = 0 ; i < (C_ij_bis_temp_0).size() ; ++i){
					long double C_1 = (C_ij_bis_temp_0)[i]/C;
					SUPERMATRICE.push_back(C_1 * G_M_bis_temp_0 );
				}
				
			}
			
			else{
				std::vector<long double> C_ij_bis_temp (C_ij(TEMP,k,ROT_CIN)); //uno per coppia
#pragma omp parallel for schedule( guided )
				for(short i = 0 ; i < (C_ij_bis_temp).size() ; ++i){
					long double C_1 = (C_ij_bis_temp)[i]/C;
					if( std::abs(C_1) > errore_moltiplicazione_matrici &&  G_M_bis_temp_0.nnz()){
						SUPERMATRICE[i].plus_assign( C_1 * (G_M_bis_temp_0) );
					}
				}
			}
			k++;
		}
		std::cerr << std::endl;
	}
	else costruisci_file = true;
	f.close();
	
	if(costruisci_file) {
		std::string Nome_dir = "matrici_G";
		int crea_dir = mkdir(Nome_dir.c_str(), 0777);
		//		if (crea_dir) std::cerr << "Directory delle matrici G esistente" << std::endl;
		//		else std::cerr << "Directory delle matrici G creata" << std::endl;
		
		f.open(s.c_str(),std::ios::out);
		std::vector<long double> C_ij_bis_temp_0 (C_ij(TEMP,0,ROT_CIN)); //uno per coppia
		MATRICI_SPARSE G_M_bis_temp_0 (G_M(TEMP,0,G,LISTA,K));
		f << boost::numeric::ublas::io::sparse(G_M_bis_temp_0) << std::endl;
		for(short i = 0 ; i < (C_ij_bis_temp_0).size() ; ++i){
			long double C_1 = (C_ij_bis_temp_0)[i]/C;
			SUPERMATRICE.push_back(C_1 * G_M_bis_temp_0 );
		}
		
		
		//faccio la somma sui [K] -> potrei lasciare tutto esteso?
		for(short k = 1 ; k < TEMP.D() ; ++k){
			std::vector<long double> C_ij_bis_temp (C_ij(TEMP,k,ROT_CIN)); //uno per coppia
			MATRICI_SPARSE G_M_bis_temp (G_M(TEMP,k,G,LISTA,K));
			f << boost::numeric::ublas::io::sparse(G_M_bis_temp) << std::endl;
#pragma omp parallel for schedule( guided )
			for(short i = 0 ; i < (C_ij_bis_temp).size() ; ++i){
				long double C_1 = (C_ij_bis_temp)[i]/C;
				if( std::abs(C_1) > errore_moltiplicazione_matrici &&  G_M_bis_temp.nnz()){
					SUPERMATRICE[i].plus_assign( C_1 * (G_M_bis_temp) );
				}
			}
		}		
		
		f.close();
	}
	
	
	return SUPERMATRICE;
}

inline V_MATRICI_SPARSE G_M_ij(const unsigned short & n,
							   const Accoppiamenti_Iperangolari & TEMP,
							   const HyperInt & G,
							   Integrazioni_3_P * LISTA,
							   const Accoppiamenti_Iperangolari_K_max & K,
							   const std::vector< vettore_ldouble > & ROT_CIN,
							   const unsigned short & dimensione_calcolata){
	unsigned short N = TEMP.N();
	V_MATRICI_SPARSE SUPERMATRICE;
	unsigned short numero_stati_aggiunti = K.D() - dimensione_calcolata;
	V_MATRICI_SPARSE SUPERMATRICE_0;
	hyper H0 (N,n); // PB
	long double C = (H0.RIPELLE( N - 1, 0. ));
	
	MATRICI_SPARSE G_M_bis_temp_0;
	
	std::fstream f;
	f.precision(std::numeric_limits<long double>::digits10);
	//	std::string s = "./matrici_G/G_M_base_" + K.fileD() + "_Kg_" + scrivi(TEMP.K())+ "_d_" + scrivi(dimensione_calcolata) + " .txt";
	std::string s0 = "./matrici_G/G_M_base_" + K.fileD(dimensione_calcolata) + "_Kg_" + scrivi(TEMP.K()) + ".txt";
	std::string s = "./matrici_G/G_M_base_" + K.fileD() + "_Kg_" + scrivi(TEMP.K()) + ".txt";
	
	//vedo se esiste
	bool costruisci_file=false;
	f.open(s.c_str(),std::ios::in);
	if (f) {
		unsigned long k=0;
		while (f >> boost::numeric::ublas::io::sparse(G_M_bis_temp_0)){
			if (k==0) {
				std::vector<long double> C_ij_bis_temp_0 (C_ij(TEMP,0,ROT_CIN)); //uno per coppia
				for(short i = 0 ; i < (C_ij_bis_temp_0).size() ; ++i){
					long double C_1 = (C_ij_bis_temp_0)[i]/C;
					SUPERMATRICE.push_back(C_1 * G_M_bis_temp_0 );
				}
				
			}
			
			else{
				std::vector<long double> C_ij_bis_temp (C_ij(TEMP,k,ROT_CIN)); //uno per coppia
#pragma omp parallel for schedule( guided )
				for(short i = 0 ; i < (C_ij_bis_temp).size() ; ++i){
					long double C_1 = (C_ij_bis_temp)[i]/C;
					if( std::abs(C_1) > errore_moltiplicazione_matrici &&  G_M_bis_temp_0.nnz()){
						SUPERMATRICE[i].plus_assign( C_1 * (G_M_bis_temp_0) );
					}
				}
			}
			k++;
		}
	}
	else costruisci_file = true;
	f.close();
	
	if(costruisci_file) {
		//leggi da un file le matrici da aggiornare SUPERMATRICE_0
		f.open(s0.c_str(),std::ios::in);
		if (f) {
			unsigned long k=0;
			while (f >> boost::numeric::ublas::io::sparse(G_M_bis_temp_0)){
				SUPERMATRICE_0.push_back(G_M_bis_temp_0);
				k++;
			}
			
		}
		else{
			std::cerr<<"**********************************" << std::endl;
			std::cerr<<"Grave Errore!! " << std::endl;
			std::cerr<<"**********************************" << std::endl;
		}
		f.close();
		//estendi la supermatrice_0
		for(int isupm = 0; isupm < SUPERMATRICE_0.size(); ++isupm){
			SUPERMATRICE_0[isupm] = accresci(SUPERMATRICE_0[isupm],numero_stati_aggiunti);
		}
		//ho esteso; calcolo le aggiunte e scrivo il file aggiornato
		f.open(s.c_str(),std::ios::out);
		std::vector<long double> C_ij_bis_temp_0 (C_ij(TEMP,0,ROT_CIN)); //uno per coppia
		MATRICI_SPARSE G_M_bis_temp_0 (G_M(TEMP,0,G,LISTA,K,dimensione_calcolata));
		SUPERMATRICE_0[0].plus_assign(G_M_bis_temp_0);
		f << boost::numeric::ublas::io::sparse(SUPERMATRICE_0[0]) << std::endl;
		for(short i = 0 ; i < (C_ij_bis_temp_0).size() ; ++i){
			long double C_1 = (C_ij_bis_temp_0)[i]/C;
			SUPERMATRICE.push_back(C_1 * SUPERMATRICE_0[0] );
		}
		
		
		//faccio la somma sui [K] -> potrei lasciare tutto esteso?
		for(short k = 1 ; k < TEMP.D() ; ++k){
			std::vector<long double> C_ij_bis_temp (C_ij(TEMP,k,ROT_CIN)); //uno per coppia
			MATRICI_SPARSE G_M_bis_temp (G_M(TEMP,k,G,LISTA,K,dimensione_calcolata));
			SUPERMATRICE_0[k].plus_assign(G_M_bis_temp);
			f << boost::numeric::ublas::io::sparse(SUPERMATRICE_0[k]) << std::endl;
#pragma omp parallel for schedule( guided )
			for(short i = 0 ; i < (C_ij_bis_temp).size() ; ++i){
				long double C_1 = (C_ij_bis_temp)[i]/C;
				if( std::abs(C_1) > errore_moltiplicazione_matrici &&  SUPERMATRICE_0[k].nnz()){
					SUPERMATRICE[i].plus_assign( C_1 * SUPERMATRICE_0[k] );
				}
			}
		}		
		
		f.close();
	}
	
	return SUPERMATRICE;
}

#endif