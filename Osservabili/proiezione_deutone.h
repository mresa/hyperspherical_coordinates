/*
 *  proiezione_deutone.h
 *  Tesi
 *
 *  Created by Marco Resa on 25/05/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef PR_D_H
#define PR_D_H
#include "tipi_definiti.h"
#include "algebra_matrici.h"
#include <boost/numeric/ublas/io.hpp>
#include <vector>
#include "accoppiamenti_iperangolari.h"
#include "jacobi.h"
#include "hyper.h"


//dato lo stato di 3 corpi, sono utili solo gli indici che hanno l_1 e l_2 = 0: ritorna gli indici dello sviluppo utili
std::vector<unsigned long> indici_utili_K_deutone(const Accoppiamenti_Iperangolari_K_max & KK){
	std::vector<unsigned long> indici;
	for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
		if(KK.l(k_i)[0] == 0 && KK.l(k_i)[1] == 0){
			indici.push_back(k_i);
		}
	}
	cerr_vec(indici);
	std::cerr << std::endl;
	return indici;
}

//matrice degli indici utili dati gli stati di 3 corpi
MATRICI A_mn(const Accoppiamenti_Iperangolari_K_max & KK,
			 const unsigned short & dimensione_m,
			 const VETTORI & mygen){
	std::vector<unsigned long> indici_k(indici_utili_K_deutone(KK));
	MATRICI A_mn_temp(dimensione_m,indici_k.size());
	for (unsigned long m = 0 ; m<A_mn_temp.size1(); ++m) {
		for (unsigned long n = 0 ; n<A_mn_temp.size2(); ++n) {
			std::cerr << "n = " << n << " m = "<< m << std::endl;
			std::cerr << "indici_k[n] = " << indici_k[n] << std::endl;
			std::cerr << indici_k[n]*m << std::endl;
			std::cerr << mygen(A_mn_temp.size2()*indici_k[n] + m) << std::endl;

			A_mn_temp(m,n)=mygen(A_mn_temp.size2()*indici_k[n] + m);
		}
	}
	return A_mn_temp;
}


long double integrando_deutone(const long double & beta, 
							   const unsigned short & m, 
							   const unsigned short & n, 
							   const unsigned short & m_d, 
							   const long double & Zeta, 
							   const long double & Ypsilon){
	long double yps = beta*Ypsilon*beta*Ypsilon/4.;
	long double zed = Zeta*Zeta;
	
	long double risultato = (
							 pow(Zeta, 2)	* 
							 laguerre ( m_d , 2 , 2.*Zeta )	*
							 laguerre ( m , 5 , 2.*sqrt(zed + yps))	*
							 exp ( -sqrt(zed + yps))	*
							 jacobi( n , 0.5 , 0.5 , (zed - yps)/(zed + yps) )
							 );
	if (risultato != risultato) {
		risultato = 0.;
		std::cerr << "errore in integrando_deutone" << std::endl;
	}
	return risultato;
}

long double costanti_integrando_deutone(const long double & beta, 
										const unsigned short & m, 
										const unsigned short & n, 
										const unsigned short & m_d){
	
	long double risultato = 8.*sqrt(
									boost::math::tgamma_delta_ratio(m + 1, 5) * 
									boost::math::tgamma_delta_ratio(m_d + 1, 2) * 
									pow(beta, 3)
									)*tre_corpi_normalizzazione(n,0,0);
	return risultato;
}

long double integrale_z_deutone(const Gauss_Laguerre & G, 
								const long double & beta, 
								const unsigned short & m, 
								const unsigned short & n, 
								const unsigned short & m_d, 
								const long double & Ypsilon){
	long double risultato = 0.;
#pragma omp parallel for reduction(+:risultato) schedule( guided )
	for (int k = 0; k < G.N() ; ++k){
		long double incremento = G.W(k) * integrando_deutone(beta,m,n,m_d,G.X(k),Ypsilon);
		if (incremento != incremento){
			incremento = 0.;	
			std::cerr << "errore in integrale_z_deutone" << std::endl;
		} 
		risultato += incremento;
	}
	risultato*=costanti_integrando_deutone(beta,m,n,m_d);
	if (risultato != risultato){
		risultato = 0.;	
		std::cerr << "errore in integrale_z_deutone - tot" << std::endl;
	} 	
	//std::cerr << "integrale_z_deutone = " << risultato << std::endl;
	return risultato;
}

long double integrale_deutone_y(const Gauss_Laguerre & G, 
							   const long double & beta, 
							   const VETTORI & c_m_d, 
							   const MATRICI & A_mn, 
							   const long double & Ypsilon){
	long double risultato = 0.;
#pragma omp parallel for reduction(+:risultato) schedule( guided )
	//	for (int k = 0; k < G.N() ; ++k){
	//	long double incremento = 0.;
	for (int m = 0; m < A_mn.size1() ; ++m)
		for (int n = 0; n < A_mn.size2() ; ++n){
			long double incremento_mn = 0.;
			for (int m_d = 0; m_d < c_m_d.size() ; ++m_d){
				incremento_mn += c_m_d(m_d)*integrale_z_deutone(G,beta,m,n,m_d,Ypsilon);
			}
			if (incremento_mn != incremento_mn){
				incremento_mn = 0.;	
				std::cerr << "errore in integrale_deutone_y" << std::endl;
			}
			else incremento_mn*=A_mn(m,n);
			if (incremento_mn != incremento_mn){
				incremento_mn = 0.;	
				std::cerr << "errore in integrale_deutone_y" << std::endl;
			}
			//incremento += incremento_mn;
			//	}
			risultato += incremento_mn;
		}
	return risultato;
}


#endif
