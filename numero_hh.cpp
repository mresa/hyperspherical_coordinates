/*
 *  numero_hh.cpp
 *  Tesi
 *
 *  Created by Marco Resa on 25/04/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "preambolo.h"
#include "gauss-laguerre.h"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include "algebra_matrici.h"
#include <boost/numeric/ublas/io.hpp>
#include <assert.h>
#include <ietl/traits.h>
#include "accoppiamenti_iperangolari.h"
#include "interfaces_ublas.h"
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>


void init_main_0(const unsigned short & A,
				 const unsigned short & K_min,
				 const unsigned short & K_max,
				 const unsigned short & P,
				 const unsigned short modo_clust,
				 const unsigned short & P_orb){
	
	unsigned short N = A - 1;
	std::string s, s_1;
	std::string s_k,s_cluster;
	std::string s_ext_txt = ".txt";
	std::string s_ext_dat = ".dat";
	std::string s_k_max = "_K_max_" + scrivi(K_max);
	
	std::vector <unsigned short> clus;
	
	if (modo_clust && A==4) {
		clus.push_back(2);
		clus.push_back(2);
		s_cluster = "_Modo_H";
	}
	else {
		clus.push_back(A);
		s_cluster = "";
	}

	Accoppiamenti_Iperangolari_K_max KK (N, 0, 0, K_min,P,clus,P_orb);
	s_k = KK.filename();
	{
		s = "numero_hh_";
		s += s_k + s_k_max + s_cluster  + s_ext_dat;
		s_1 = "numero_hh_";
		s_1 += s_k + s_k_max + s_cluster  + s_ext_txt;
	}
	
		
	std::ofstream myfile;
	std::ofstream myfile_1;
	{
		myfile.open (s.c_str());
		myfile_1.open (s_1.c_str());
		myfile.precision(std::numeric_limits<long double>::digits10); 		
		myfile_1.precision(std::numeric_limits<long double>::digits10); 		
		std::cerr.precision(10);
		std::cout.precision(10);
	}
	
	{		
		myfile_1 << "Numero di HH per " << A << " nucleoni ";
		myfile_1 << std::endl;
	}
	
	boost::timer t6;
	unsigned short incremento = 1;

	for (unsigned short K_max_temp = K_min; K_max_temp <= K_max; K_max_temp+=incremento) {
				
		boost::timer t;
		
		{	
			myfile_1 << std::endl;
			myfile_1 << std::endl;
			myfile_1 << "**************************************************************************";
			myfile_1 << std::endl;
			myfile_1 << " K_max_temp = " << K_max_temp << std::endl;
			std::cerr << " K_max_temp = " << K_max_temp << std::endl;
		}
		KK.PLUS(K_max_temp - KK.K_M(),P,clus);
		if (KK.D()>0) {
			{
				myfile_1 << " N_hh = " << KK.D() << std::endl;
				std::cerr << " N_hh = " << KK.D() << std::endl;
				std::cerr << "Accoppiamenti_Iperangolari_K_max  T = " << t.elapsed() << std::endl;
				myfile << "Accoppiamenti_Iperangolari_K_max  T = " << t.elapsed() << std::endl;
				for (unsigned long i = 0; i < KK.D() ; ++i){
					myfile << i+1 << "	K = " << KK.k(i) << "	 n = [ ";
					for(unsigned short j = 0; j < KK.n(i).size() ; ++j){
						myfile << KK.n(i)[j];
						myfile << " ";
					}
					myfile << "]   l = [ ";
					for(unsigned short j = 0; j < KK.l(i).size() ; ++j){
						myfile << KK.l(i)[j];
						myfile << " ";
					}
					
					myfile << "]  L = [ ";
					for(unsigned short j = 0; j < KK.L(i).size() ; ++j){
						myfile << KK.L(i)[j];
						myfile << " ";
					}
					
					myfile << "]";
					
					myfile << "   parita = [ ";
					for(unsigned short j = 0; j < KK.P_l(i).size() ; ++j){
						myfile << KK.P_l(i)[j];
						myfile << " ";
					}
					myfile << "]	 somma quadra dei coefficienti dello sviluppo = ";
					myfile << quad_sum(KK.cm(i));				
					myfile << std::endl;
				}
				myfile << "Dimensione della base = " << KK.D() << " combinazioni." << std::endl;
			}
		}
	}
	std::cerr << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile.close();
	myfile_1.close();
	
}

int main (int argc, char *argv[]) {
	
	MKLVersion ver;
	
	// Print information on CPU optimization in effect
	
	MKLGetVersion(&ver);
	
	printf("Processor optimization: %s\n",ver.Processor);
	omp_set_num_threads(8);
	mkl_set_dynamic(true);
	
	unsigned short A = ((unsigned short) std::strtoul(argv[1],(char **)0, 10));
	unsigned short K_min = ((unsigned short) std::strtoul(argv[2],(char **)0, 10));
	unsigned short K_max = ((unsigned short) std::strtoul(argv[3],(char **)0, 10));
	unsigned short P = ((unsigned short) std::strtoul(argv[4],(char **)0, 10));
	unsigned short Parita_Orbitale = ((unsigned short) std::strtoul(argv[5],(char **)0, 10));; //0 = pari ; 1 = dispari ; 2 = tutte
	unsigned short modo_cluster = 0;
	if (A==4) {
		if (argc == 7) {
			modo_cluster = ((unsigned short) std::strtoul(argv[6],(char **)0, 10)); //mettere a 1 per avere modo H sui 4 corpi
		}
	}
	init_main_0(A,K_min,K_max,P,modo_cluster,Parita_Orbitale);
	return true;
}



