/*
 *  energia_3_corpi_orbitale
 *  Tesi
 *
 *  Created by Marco Resa on 11/02/10.
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
#include "rotazione_iperangolare.h"
#include "accoppiamenti_parita.h"
#include "init_corpi.h"
#include "energia.h"
#include "interfaces_ublas.h"
#include <ietl/vectorspace.h>
#include <ietl/lanczos.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>
#include "autostati.h"


void init_main_0(const unsigned short & A,
				 const unsigned short & m_max,
				 const unsigned short & K_min,
				 const unsigned short & K_max,
				 const Init_Corpi_L & Corpi,
				 const unsigned short & P,
				 const bool & modo_clust,
				 const unsigned short & P_orb){
	
	unsigned short ordine = 200;
	unsigned short N = A - 1;
	
	std::string s, s_0, s_1 , s_l, s_misurabili;
	std::string s_0_1, s_misurabili_1;
	std::string s_k,s_cluster,s_masse;
	std::string s_ext_txt = "_ecc.txt";
	std::string s_ext_dat = "_ecc.dat";
	std::string s_k_max = "_K_max_" + scrivi(K_max);
	std::string s_corpi = "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max);
	std::string s_ecc = "_ecc1";

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
	if(Corpi.M_u()){
		s_masse = "_M_u";
	}
	else {
		s_masse = "";
	}
	
	Accoppiamenti_Iperangolari_K_max KK (N, 0, 0, K_min,P,clus,P_orb);
	s_k = KK.filename();
	{
		s = "Nuova_E_";
		s_0 = "Nuova_Tabella_E_";
		s_1 = "Nuova_E_";
		s_l = "Livelli_E_";
		s_misurabili = "Misurabili_";
		s_0_1 = "Nuova_Tabella_E_";
		s_misurabili_1 = "Misurabili_";
	}
	{
		s += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ext_dat;
		s_0 += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ext_txt;
		s_1 += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ext_txt;
		s_l += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ext_txt;
		s_misurabili += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ext_txt;
		s_0_1 += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ecc + s_ext_txt;
		s_misurabili_1 += s_corpi + s_k + s_k_max + s_cluster + s_masse + s_ecc + s_ext_txt;
	}
	std::string s_vec,s_vec_temp,s_vec_temp_c;
	{
		s_vec = "./lista_vettori/Vec_ecc_";
		s_vec += s_corpi + s_k + s_k_max + s_cluster + s_masse;
	}
	std::fstream myfile_vec;
	myfile_vec.precision(10); 		
	
	std::ofstream myfile;
	std::ofstream myfile_1;
	std::ofstream myfile_tabella;
	std::ofstream myfile_livelli;
	std::ofstream myfile_misurabili;
	std::ofstream myfile_tabella_1;
	std::ofstream myfile_misurabili_1;
	{
		myfile.open (s.c_str());
		myfile_1.open (s_1.c_str());
		myfile_tabella.open (s_0.c_str());
		myfile_tabella_1.open (s_0_1.c_str());
		myfile_livelli.open (s_l.c_str());
		myfile_misurabili.open (s_misurabili.c_str());
		myfile_misurabili_1.open (s_misurabili_1.c_str());
		myfile.precision(std::numeric_limits<long double>::digits10); 		
		myfile_1.precision(std::numeric_limits<long double>::digits10); 		
		myfile_tabella.precision(10); 		
		myfile_tabella_1.precision(10); 		
		myfile_livelli.precision(10); 		
		myfile_misurabili.precision(10); 		
		myfile_misurabili_1.precision(10); 		
		std::cerr.precision(10);
		std::cout.precision(10);
	}
	
	{		
		myfile_1 << "Binding Energy per " << A << " nucleoni : ";
		myfile_tabella << "Binding Energy per " << A << " nucleoni : ";
		myfile_tabella_1 << "Binding Energy per il primo eccitato per " << A << " nucleoni : ";
		myfile_livelli << "Primi Livelli Energetici per " << A << " nucleoni : ";
		for (short i = 0; i<A; ++i) {
			myfile_1<< Corpi.CORPI(i)<< " ";
			myfile_tabella<< Corpi.CORPI(i)<< " ";
			myfile_tabella_1<< Corpi.CORPI(i)<< " ";
			myfile_livelli<< Corpi.CORPI(i)<< " ";
		}
		myfile_1 << std::endl;
		myfile_tabella << std::endl;
		myfile_tabella_1 << std::endl;
		myfile_livelli << std::endl;
		myfile_1 << " m_max = " << m_max << " K_max = " << K_max << "	P = " << P << std::endl;
		myfile_tabella << " m_max = " << m_max << " K_max = " << K_max << "	P = " << P << std::endl;
		myfile_tabella_1 << " m_max = " << m_max << " K_max = " << K_max << "	P = " << P << std::endl;
		myfile_livelli << " m_max = " << m_max << " K_max = " << K_max << "	P = " << P << std::endl;
	}
	if(Corpi.M_u() && A==3)	{
		myfile_tabella << "K_max	N_hh	Parita'		E		E_ppp		E_ppn		E_nnp		E_nnn		E_c" << std::endl;
		myfile_tabella_1 << "K_max	N_hh	Parita'		E		E_ppp		E_ppn		E_nnp		E_nnn		E_c" << std::endl;
	}
	else{
		myfile_tabella << "K_max	N_hh	Parita'		E		E_c" << std::endl;
		myfile_tabella_1 << "K_max	N_hh	Parita'		E		E_c" << std::endl;
	}
	boost::timer t6;
	unsigned short incremento = 2;
	Gauss_Laguerre G(ordine);
	HyperInt HI(ordine/2);
	Energia_radiale ER ( N, G, HI, K_min, m_max, Corpi.C(), Corpi.M(), Corpi.M_R(), clus);
	Integrazioni_3_P INT_3P(N);
	Energia EN ( HI, &INT_3P , ER, KK );
	
	
	//lancsoz
	VETTORI mygen(m_max + 1);
	VETTORI mygen_c(m_max + 1);
	typedef boost::lagged_fibonacci607 Gen_i;
	typedef ietl::vectorspace<VETTORI> Vecspace;
	Gen_i mygen_i;
	bool primopasso = true;
	tipo_matrici rel_tol = 500000000*std::numeric_limits<double>::epsilon();
	tipo_matrici abs_tol = 1000000*std::pow(std::numeric_limits<double>::epsilon(),2./3);
	int n_lowest_eigenval = 1;	
	//lancsoz
	
	
	for (unsigned short K_max_temp = K_min; K_max_temp <= K_max; K_max_temp+=incremento) {
		
		if (K_max_temp == 20) incremento = 10;
		
		boost::timer t;
		
		{	
			myfile_1 << std::endl;
			myfile_1 << std::endl;
			myfile_1 << "**************************************************************************";
			myfile_1 << std::endl;
			myfile_1 << " K_max_temp = " << K_max_temp << std::endl;
			std::cerr << " K_max_temp = " << K_max_temp << std::endl;
		}
		{	
			myfile_livelli << std::endl;
			myfile_livelli << std::endl;
			myfile_livelli << "**************************************************************************";
			myfile_livelli << std::endl;
			myfile_livelli << " K_max_temp = " << K_max_temp << std::endl;
		}	
		KK.PLUS(K_max_temp - KK.K_M(),P,clus);
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
				
				myfile << "]    parita = [ ";
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
		
		if (KK.D()>0) {
			ER.PLUS((K_max_temp - ER.N_M()),G,HI);
			boost::timer t1;
			EN.PLUS(HI, &INT_3P, ER, KK);
			
			{
				
				std::cerr << "Energia  T = " << t1.elapsed() << std::endl;
				myfile << "Energia  T = " << t1.elapsed() << std::endl;
				
			}
			Accoppiamenti_Parita BASE(KK,ER.D());
			int max_iter = 5*BASE.D();  
			
			std::cerr << "BASE.D() "<< BASE.D() << std::endl;
			
			/*if(BASE.D()>11500){
				
				std::cerr << " Sto popolando più di 2 Gb, mi fermo" << std::endl;
				
				if (primopasso) {
					mygen.resize(BASE.D());
					mygen_c.resize(BASE.D());
					ietl::generate_begin(mygen,mygen_i);
					ietl::generate_begin(mygen_c,mygen_i);
					primopasso = false;
				}
				else{
					ietl::regenerate_begin(mygen,mygen_i,BASE.D());					
					ietl::regenerate_begin(mygen_c,mygen_i,BASE.D());					
				}
				
				myfile_tabella << K_max_temp << "	";
				myfile_tabella << KK.D() << "	";
				myfile_tabella << scrivi(P) <<  "	";
				myfile_misurabili << std::endl << "K_max_temp " << K_max_temp << std::endl;
				myfile_misurabili << std::endl;
				myfile_misurabili << "P = "<< scrivi(P) << std::endl;
				
				{
					boost::timer t3;
					Vecspace vec(BASE.D());
					EN.N_COMPUTA_COULOMB();
					mkl_free_buffers();
					ietl::lanczos<Energia,Vecspace> lanczos (EN,vec);
					// Creation of an iteration object:  
					std::cerr << "Calcolo iterativo del primo autovalore - Autovalori E" << std::endl;
					std::cerr << "-----------------------------------"<< std::endl;
					myfile_1 << "Calcolo iterativo del primo autovalore - Autovalori E" << std::endl;
					myfile_1 << "-----------------------------------"<< std::endl;
					std::vector<tipo_matrici> eigen;
					std::vector<tipo_matrici> err;
					std::vector<int> multiplicity;  
					
					
					
					std::cerr << "     " ;
					
					ietl::lanczos_iteration_nlowest<tipo_matrici> iter (max_iter, n_lowest_eigenval, rel_tol, abs_tol);
					try{
						lanczos.calculate_eigenvalues(iter,mygen);
						eigen = lanczos.eigenvalues();
						err = lanczos.errors();
						multiplicity = lanczos.multiplicities();
						std::cerr<< std::endl << "numero di iterazioni "<<iter.iterations()<< std::endl;
						myfile_1<<"numero di iterazioni "<<iter.iterations()<< std::endl;
					}
					catch (std::runtime_error& e) {
						std::cerr << e.what() << "\n";
					}	
					mkl_free_buffers();
					
					{
						myfile_1 << std::endl;
						myfile_1 <<  "Autovalori E" << std::endl;
						std::cout << "#        autovalore            errore         molteplicità\n";  
						myfile_1 << "#        autovalore            errore         molteplicità\n";  
						for (int i=0;(i<eigen.size())&&(i<n_lowest_eigenval+10);++i) {
							myfile_1 << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";
						}
						for (int i=0;(i<eigen.size())&&(i<n_lowest_eigenval+3);++i) {
							std::cerr << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";
						}
						myfile_1 << std::endl;
						
						myfile_tabella << eigen[0] << "	" ;
						{
							myfile_livelli <<  std::endl << "Autovalori E"<<std::endl;
							for (int i=0; i< eigen.size() ||  i < 10 ; ++i)	{
								myfile_livelli <<  eigen[i]<< "		" ;
							}
							myfile_livelli << std::endl;
						}						
						
					}
					
					// call of eigenvectors function follows:   
					{
						boost::timer t4;
						std::cerr << "\nRicerca degli autovettori per il primo autovalore:\n\n";  
						myfile_1 << "\nRicerca degli autovettori per il primo autovalore:\n\n";  
						std::vector<tipo_matrici>::iterator start = eigen.begin();  
						std::vector<tipo_matrici>::iterator end = eigen.begin()+1;
						std::vector<VETTORI> eigenvectors;
						ietl::Info<tipo_matrici> info; // (m1, m2, ma, eigenvalue, residualm, status).
						std::cerr << "     " ;
						
						try {
							lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen,max_iter); 
							std::cerr<< std::endl << "numero di iterazioni "<<iter.iterations()<< std::endl;
							myfile_1<<"numero di iterazioni "<<iter.iterations()<< std::endl;
						}
						catch (std::runtime_error& e) {
							std::cout << e.what() << "\n";
						}  
						mkl_free_buffers();
						
						std::cout << "Autovettori calcolati:\n\n"; 
						
						
						std::cout << " Informazioni sugli autovettori:\n\n";
						for(int i = 0; i < info.size(); i++) {
							std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
							<< info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
							<< i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
							<< info.residual(i) << " error_info(" << i+1 << "): "
							<< info.error_info(i) << "\n\n"<< " normalizzazione = " << boost::numeric::ublas::norm_2 (eigenvectors[0]) << std::endl;
							myfile_1 << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
							<< info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
							<< i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
							<< info.residual(i) << " error_info(" << i+1 << "): "
							<< info.error_info(i) << "\n\n" << " normalizzazione = " << boost::numeric::ublas::norm_2 (eigenvectors[0]) << std::endl;
						}
						
						std::cerr << "Autovalori indiretta  T = " << t4.elapsed() << std::endl;
						myfile << "Autovalori indiretta  T = " << t4.elapsed() << std::endl;
						
						if (info.error_info(0)==0){
							mygen.assign(eigenvectors[0]);
							
							VETTORI autostato(eigenvectors[0]);
							
							boost::timer t5;
							std::cerr << " Calcolo dei misurabili per l'autovettore:\n\n";
							myfile_1 << " Calcolo dei misurabili per l'autovettore:\n\n";
							{
								
								VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
								MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
								tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
								tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
								
								std::cout << "r_i^2 =  " << v1 << std::endl;
								std::cout << "r_ij^2 =  " << v2 << std::endl;
								std::cout << std::endl;
								
								std::cout << "T =  " << v3 << std::endl;
								std::cout << "V =  " << v4 << std::endl;
								std::cout << std::endl;
								
								myfile_1 << "r_i^2 =  " << v1 << std::endl;
								myfile_1 << "r_ij^2 =  " << v2 << std::endl;
								myfile_1 << "T =  " << v3 << std::endl;
								myfile_1 << "V =  " << v4 << std::endl;
								myfile_1 << std::endl;
								
								myfile_misurabili << "r_i^2 =  " << v1 << std::endl;
								myfile_misurabili << "r_ij^2 =  " << v2 << std::endl;
								myfile_misurabili << "T =  " << v3 << std::endl;
								myfile_misurabili << "V =  " << v4 << std::endl;
								myfile_misurabili << std::endl;
								
								if(Corpi.M_u() && A==3)	{
									tipo_matrici v5 = valor_medio(EN.MATRICI_NABLA(),autostato);
									
									myfile_tabella << eigen[0] + (m_A + m_B)*v3 << "	";
									myfile_tabella << eigen[0] + (m_A + m_B)*v3 + (4./3.)*m_B*v5<< "	";
									myfile_tabella << eigen[0] + (m_A - m_B)*v3 - (4./3.)*m_B*v5<< "	";
									myfile_tabella << eigen[0] + (m_A - m_B)*v3 << "	";
								}
								
							}
							mkl_free_buffers();
							
							std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
							myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
						}
					}
					
					
					std::cerr << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					myfile << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					
				}
				{
					boost::timer t3;
					Vecspace vec(BASE.D());
					EN.COMPUTA_COULOMB();
					mkl_free_buffers();
					ietl::lanczos<Energia,Vecspace> lanczos (EN,vec);
					// Creation of an iteration object:  
					std::cerr << "Calcolo iterativo del primo autovalore - Autovalori E_C" << std::endl;
					std::cerr << "-----------------------------------"<< std::endl;
					myfile_1 << "Calcolo iterativo del primo autovalore - Autovalori E_C" << std::endl;
					myfile_1 << "-----------------------------------"<< std::endl;
					std::vector<tipo_matrici> eigen;
					std::vector<tipo_matrici> err;
					std::vector<int> multiplicity;  
					
					
					
					std::cerr << "     " ;
					
					ietl::lanczos_iteration_nlowest<tipo_matrici> iter (max_iter, n_lowest_eigenval, rel_tol, abs_tol);
					try{
						lanczos.calculate_eigenvalues(iter,mygen_c);
						eigen = lanczos.eigenvalues();
						err = lanczos.errors();
						multiplicity = lanczos.multiplicities();
						std::cerr<< std::endl << "numero di iterazioni "<<iter.iterations()<< std::endl;
						myfile_1<<"numero di iterazioni "<<iter.iterations()<< std::endl;
					}
					catch (std::runtime_error& e) {
						std::cerr << e.what() << "\n";
					}	
					mkl_free_buffers();
					
					{
						myfile_1 << std::endl;
						myfile_1 <<  "Autovalori E_C" << std::endl;
						std::cout << "#        autovalore            errore         molteplicità\n";  
						myfile_1 << "#        autovalore            errore         molteplicità\n";  
						for (int i=0;(i<eigen.size())&&(i<n_lowest_eigenval+10);++i) {
							myfile_1 << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";
						}
						for (int i=0;(i<eigen.size())&&(i<n_lowest_eigenval+3);++i) {
							std::cerr << i << "\t" << eigen[i] << "\t" << err[i] << "\t" << multiplicity[i] << "\n";
						}
						myfile_1 << std::endl;
						
						myfile_tabella << eigen[0] << std::endl;
						{
							myfile_livelli <<  std::endl << "Autovalori E_c"<<std::endl;
							for (int i=0; i< eigen.size() ||  i < 10 ; ++i)	{
								myfile_livelli <<  eigen[i]<< "		" ;
							}
							myfile_livelli << std::endl;
						}						
						
					}
					
					// call of eigenvectors function follows:   
					{
						boost::timer t4;
						std::cerr << "\nRicerca degli autovettori per il primo autovalore:\n\n";  
						myfile_1 << "\nRicerca degli autovettori per il primo autovalore:\n\n";  
						std::vector<tipo_matrici>::iterator start = eigen.begin();  
						std::vector<tipo_matrici>::iterator end = eigen.begin()+1;
						std::vector<VETTORI> eigenvectors;
						ietl::Info<tipo_matrici> info; // (m1, m2, ma, eigenvalue, residualm, status).
						std::cerr << "     " ;
						
						try {
							lanczos.eigenvectors(start,end,std::back_inserter(eigenvectors),info,mygen_c,max_iter); 
							std::cerr<< std::endl << "numero di iterazioni "<<iter.iterations()<< std::endl;
							myfile_1<<"numero di iterazioni "<<iter.iterations()<< std::endl;
						}
						catch (std::runtime_error& e) {
							std::cout << e.what() << "\n";
						}  
						mkl_free_buffers();
						
						std::cout << "Autovettori calcolati:\n\n"; 
						
						
						std::cout << " Informazioni sugli autovettori:\n\n";
						for(int i = 0; i < info.size(); i++) {
							std::cout << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
							<< info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
							<< i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
							<< info.residual(i) << " error_info(" << i+1 << "): "
							<< info.error_info(i) << "\n\n"<< " normalizzazione = " << boost::numeric::ublas::norm_2 (eigenvectors[0]) << std::endl;
							myfile_1 << " m1(" << i+1 << "): " << info.m1(i) << ", m2(" << i+1 << "): "
							<< info.m2(i) << ", ma(" << i+1 << "): " << info.ma(i) << " eigenvalue("
							<< i+1 << "): " << info.eigenvalue(i) << " residual(" << i+1 << "): "
							<< info.residual(i) << " error_info(" << i+1 << "): "
							<< info.error_info(i) << "\n\n" << " normalizzazione = " << boost::numeric::ublas::norm_2 (eigenvectors[0]) << std::endl;
						}
						
						std::cerr << "Autovalori indiretta  T = " << t4.elapsed() << std::endl;
						myfile << "Autovalori indiretta  T = " << t4.elapsed() << std::endl;
						
						if (info.error_info(0)==0){
							mygen_c.assign(eigenvectors[0]);
							VETTORI autostato(eigenvectors[0]);
							
							boost::timer t5;
							std::cerr << " Calcolo dei misurabili per l'autovettore:\n\n";
							myfile_1 << " Calcolo dei misurabili per l'autovettore:\n\n";
							{
								
								VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
								MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
								tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
								tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
								tipo_matrici v5 = valor_medio(EN.MATRICI_V_C(),autostato);
								
								std::cout << "r_i^2 =  " << v1 << std::endl;
								std::cout << "r_ij^2 =  " << v2 << std::endl;
								std::cout << std::endl;
								
								std::cout << "T =  " << v3 << std::endl;
								std::cout << "V =  " << v4 << std::endl;
								std::cout << "V_c =  " << v5 << std::endl;
								std::cout << std::endl;
								
								myfile_1 << "r_i^2 =  " << v1 << std::endl;
								myfile_1 << "r_ij^2 =  " << v2 << std::endl;
								myfile_1 << "T =  " << v3 << std::endl;
								myfile_1 << "V =  " << v4 << std::endl;
								myfile_1 << "V_c =  " << v5 << std::endl;
								myfile_1 << std::endl;
								
								myfile_misurabili << "Coulomb:" << std::endl;
								myfile_misurabili << "r_i^2 =  " << v1 << std::endl;
								myfile_misurabili << "r_ij^2 =  " << v2 << std::endl;
								myfile_misurabili << "T =  " << v3 << std::endl;
								myfile_misurabili << "V =  " << v4 << std::endl;
								myfile_misurabili << "V_c =  " << v5 << std::endl;
								myfile_misurabili << "V_t =  " << (v4+v5) << std::endl;
								myfile_misurabili << std::endl;
								
							}
							mkl_free_buffers();
							std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
							myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
						}
						
					}
					std::cerr << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					myfile << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					
				}
				std::cerr << "K_MAX = " << K_max_temp << " T = " << t.elapsed() << std::endl;
				myfile << "K_MAX = " << K_max_temp << " T = " << t.elapsed() << std::endl;
				
				mkl_free_buffers();
				
				
				
			}
			
			else*/
			{
				
				{
					std::string s_K_max_temp = "_" + scrivi(K_max_temp);
					s_vec_temp = s_vec + s_K_max_temp + s_ext_txt;
				}
				{
					std::string s_K_max_temp = "_" + scrivi(K_max_temp) + "_c";
					s_vec_temp_c = s_vec + s_K_max_temp + s_ext_txt;
				}

				primopasso = false;
				mygen.resize(BASE.D());
				mygen_c.resize(BASE.D());
				
				t1.restart();
				ublas::symmetric_matrix<tipo_matrici> E_composta(composizione_matrici(EN.MATRICI_E()));
				ublas::symmetric_matrix<tipo_matrici> E_composta_c(E_composta);
				if (EN.MATRICI_V_C().size()) {
					E_composta_c+=composizione_matrici(EN.MATRICI_V_C());
				}
				std::cerr << "Matrice  T = " << t1.elapsed() << std::endl;
				myfile << "Matrice  T = " << t1.elapsed() << std::endl;
				
				ublas::matrix<tipo_matrici> autovettori(E_composta.size1(),E_composta.size2());
				ublas::vector<tipo_matrici> autovalori(E_composta.size1());
				ublas::matrix<tipo_matrici> autovettori_c(E_composta_c.size1(),E_composta.size2());
				ublas::vector<tipo_matrici> autovalori_c(E_composta_c.size1());
				unsigned short indice_stato_fondamentale = 0;
				unsigned short indice_stato_eccitato = 1;
				unsigned short indice_stato_fondamentale_c = 0;
				unsigned short indice_stato_eccitato_c = 1;
				{
					boost::timer t3;
					
					autostati(E_composta, autovettori, autovalori);
					autostati(E_composta_c, autovettori_c, autovalori_c);
					
					std::cerr << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					myfile << "Diagonalizzazione indiretta  T = " << t3.elapsed() << std::endl;
					
					myfile_tabella << K_max_temp << "	";
					myfile_tabella << KK.D() << "	";
					myfile_tabella << scrivi(P) <<  "	";
					
					myfile_tabella_1 << K_max_temp << "	";
					myfile_tabella_1 << KK.D() << "	";
					myfile_tabella_1 << scrivi(P) <<  "	";
					
					myfile_misurabili << std::endl << "K_max_temp " << K_max_temp << std::endl;
					myfile_misurabili << std::endl;
					myfile_misurabili << "P = " << scrivi(P) << std::endl;
					
					myfile_misurabili_1 << std::endl << "K_max_temp " << K_max_temp << std::endl;
					myfile_misurabili_1 << std::endl;
					myfile_misurabili_1 << "P = " << scrivi(P) << std::endl;

					boost::timer t2;
					myfile_1 << std::endl;
					myfile_1 <<  "Autovalori E" << std::endl;
					for (int i=0; i< autovalori.size() ; ++i)	{
						myfile_1 <<  autovalori(i)<< "		" ;
						if (autovalori(i) < autovalori(indice_stato_fondamentale)) indice_stato_fondamentale = i;
						else if (autovalori(i) > autovalori(indice_stato_fondamentale) && autovalori(i) < autovalori(indice_stato_eccitato)) indice_stato_eccitato = i;
						if (autovalori(i) < 0 ) std::cerr << i << "		" <<  autovalori(i) << std::endl;
					}
					myfile_1 << std::endl;
					
					myfile_tabella << autovalori(indice_stato_fondamentale) << "	" ;
					myfile_tabella_1 << autovalori(indice_stato_eccitato) << "	" ;
					
					
					
					myfile_1 << std::endl;
					myfile_1 <<  "Autovalori E coulombiano " << std::endl;
					for (int i=0; i< autovalori_c.size() ; ++i)	{
						myfile_1 <<  autovalori_c(i)<< "		" ;
						if (autovalori_c(i) < autovalori_c(indice_stato_fondamentale_c)) indice_stato_fondamentale_c = i;
						else if (autovalori_c(i) > autovalori_c(indice_stato_fondamentale_c) && autovalori_c(i) < autovalori_c(indice_stato_eccitato_c)) indice_stato_eccitato_c = i;
						if (autovalori_c(i) < 0 ) std::cerr << i << "	coulombiano	" <<  autovalori_c(i) << std::endl;
					}
					myfile_1 << std::endl;
					myfile_1 << std::endl;
					
					
					
					{
						myfile_livelli <<  std::endl << "Autovalori E"<<std::endl;
						for (int i=0; i< autovalori.size() ||  i < 16 ; ++i)	{
							myfile_livelli <<  autovalori(i)<< "		" ;
						}
						myfile_livelli << std::endl;
						myfile_livelli <<  std::endl << "Autovalori E_c"<<std::endl;
						for (int i=0; i< autovalori_c.size() ||  i < 16 ; ++i)	{
							myfile_livelli <<  autovalori_c(i)<< "		" ;
						}
						myfile_livelli << std::endl;
					}
					
					std::cerr << "Autovalori  T = " << t2.elapsed() << std::endl;
					myfile << "Autovalori  T = " << t2.elapsed() << std::endl;
					{
						VETTORI autostato(column(autovettori, indice_stato_fondamentale));
						mygen.assign(autostato);
						
						boost::timer t5;
						std::cerr << " Calcolo dei misurabili per l'autovettore:\n\n";
						myfile_1 << " Calcolo dei misurabili per l'autovettore:\n\n";
						{
							
							VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
							MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
							tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
							tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
							
							std::cout << "Autovettore " << indice_stato_fondamentale << std::endl;
							std::cout << "r_i^2 =  " << v1 << std::endl;
							std::cout << "r_ij^2 =  " << v2 << std::endl;
							std::cout << std::endl;
							
							std::cout << "T =  " << v3 << std::endl;
							std::cout << "V =  " << v4 << std::endl;
							std::cout << std::endl;
							
							myfile_1 << "Autovettore " << indice_stato_fondamentale << std::endl;
							myfile_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_1 << "T =  " << v3 << std::endl;
							myfile_1 << "V =  " << v4 << std::endl;
							myfile_1 << std::endl;
							
							myfile_misurabili << "r_i^2 =  " << v1 << std::endl;
							myfile_misurabili << "r_ij^2 =  " << v2 << std::endl;
							myfile_misurabili << "T =  " << v3 << std::endl;
							myfile_misurabili << "V =  " << v4 << std::endl;
							myfile_misurabili << std::endl;
							
							if(Corpi.M_u() && A==3)	{
								tipo_matrici v5 = valor_medio(EN.MATRICI_NABLA(),autostato);
								
								myfile_tabella << autovalori(indice_stato_fondamentale) + (m_A + m_B)*v3 << "	";
								myfile_tabella << autovalori(indice_stato_fondamentale) + (m_A + m_B)*v3 + (4./3.)*m_B*v5<< "	";
								myfile_tabella << autovalori(indice_stato_fondamentale) + (m_A - m_B)*v3 - (4./3.)*m_B*v5<< "	";
								myfile_tabella << autovalori(indice_stato_fondamentale) + (m_A - m_B)*v3 << "	";
							}
							
							mkl_free_buffers();
						}
						
						std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
						myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
					}
					{
						VETTORI autostato(column(autovettori, indice_stato_eccitato));
						myfile_vec.open(s_vec_temp.c_str(),std::ios::out);
						myfile_vec << autostato << std::endl;
						myfile_vec.close();
						
						boost::timer t5;
						std::cerr << " Calcolo dei misurabili per l'autovettore:\n\n";
						myfile_1 << " Calcolo dei misurabili per l'autovettore:\n\n";
						{
							
							VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
							MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
							tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
							tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
							
							std::cout << "Autovettore " << indice_stato_eccitato << std::endl;
							std::cout << "r_i^2 =  " << v1 << std::endl;
							std::cout << "r_ij^2 =  " << v2 << std::endl;
							std::cout << std::endl;
							
							std::cout << "T =  " << v3 << std::endl;
							std::cout << "V =  " << v4 << std::endl;
							std::cout << std::endl;
							
							myfile_1 << "Autovettore eccitato " << indice_stato_eccitato << std::endl;
							myfile_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_1 << "T =  " << v3 << std::endl;
							myfile_1 << "V =  " << v4 << std::endl;
							myfile_1 << std::endl;
							
							myfile_misurabili_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_misurabili_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_misurabili_1 << "T =  " << v3 << std::endl;
							myfile_misurabili_1 << "V =  " << v4 << std::endl;
							myfile_misurabili_1 << std::endl;
							
							if(Corpi.M_u() && A==3)	{
								tipo_matrici v5 = valor_medio(EN.MATRICI_NABLA(),autostato);
								
								myfile_tabella_1 << autovalori(indice_stato_eccitato) + (m_A + m_B)*v3 << "	";
								myfile_tabella_1 << autovalori(indice_stato_eccitato) + (m_A + m_B)*v3 + (4./3.)*m_B*v5<< "	";
								myfile_tabella_1 << autovalori(indice_stato_eccitato) + (m_A - m_B)*v3 - (4./3.)*m_B*v5<< "	";
								myfile_tabella_1 << autovalori(indice_stato_eccitato) + (m_A - m_B)*v3 << "	";
							}
							
							mkl_free_buffers();
						}
						
						std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
						myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
					}
					
					myfile_tabella << autovalori_c(indice_stato_fondamentale_c) << std::endl;
					myfile_tabella_1 << autovalori_c(indice_stato_eccitato_c) << std::endl;
					
					{
						VETTORI autostato(column(autovettori_c, indice_stato_fondamentale_c));
						mygen_c.assign(autostato);
						
						boost::timer t5;
						std::cerr << " Calcolo dei misurabili per l'autovettore - coulomb:\n\n";
						myfile_1 << std::endl;
						myfile_1 << " Calcolo dei misurabili per l'autovettore - coulomb:\n\n";
						{
							
							VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
							MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
							tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
							tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
							tipo_matrici v5 = valor_medio(EN.MATRICI_V_C(),autostato);
							
							std::cout << "Autovettore " << indice_stato_fondamentale_c << std::endl;
							std::cout << "r_i^2 =  " << v1 << std::endl;
							std::cout << "r_ij^2 =  " << v2 << std::endl;
							std::cout << std::endl;
							
							std::cout << "T =  " << v3 << std::endl;
							std::cout << "V =  " << v4 << std::endl;
							std::cout << "V_c =  " << v5 << std::endl;
							std::cout << std::endl;
							
							myfile_1 << "Autovettore " << indice_stato_fondamentale_c << std::endl;
							myfile_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_1 << "T =  " << v3 << std::endl;
							myfile_1 << "V =  " << v4 << std::endl;
							myfile_1 << "V_c =  " << v5 << std::endl;
							myfile_1 << std::endl;
							
							myfile_misurabili << "Coulomb:" << std::endl;
							myfile_misurabili << "r_i^2 =  " << v1 << std::endl;
							myfile_misurabili << "r_ij^2 =  " << v2 << std::endl;
							myfile_misurabili << "T =  " << v3 << std::endl;
							myfile_misurabili << "V =  " << v4 << std::endl;
							myfile_misurabili << "V_c =  " << v5 << std::endl;
							myfile_misurabili << "V_t =  " << (v4+v5) << std::endl;
							myfile_misurabili << std::endl;
							
							mkl_free_buffers();
						}
						
						std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
						myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
					}
					{
						VETTORI autostato(column(autovettori_c, indice_stato_eccitato_c));
						myfile_vec.open(s_vec_temp_c.c_str(),std::ios::out);
						myfile_vec << autostato << std::endl;
						myfile_vec.close();
						
						boost::timer t5;
						std::cerr << " Calcolo dei misurabili per l'autovettore - coulomb:\n\n";
						myfile_1 << std::endl;
						myfile_1 << " Calcolo dei misurabili per l'autovettore - coulomb:\n\n";
						{
							
							VETTORI v1(r_quadro(EN.MATRICI_R_QUADRO(),autostato));
							MATRICI_SIMMETRICHE_UP v2(r_delta_quadro(EN.MATRICI_R_DELTA_QUADRO(),A,autostato));
							tipo_matrici v3 = valor_medio(EN.MATRICI_T(),autostato);
							tipo_matrici v4 = valor_medio(EN.MATRICI_V(),autostato);
							tipo_matrici v5 = valor_medio(EN.MATRICI_V_C(),autostato);
							
							std::cout << "Autovettore eccitato " << indice_stato_eccitato_c << std::endl;
							std::cout << "r_i^2 =  " << v1 << std::endl;
							std::cout << "r_ij^2 =  " << v2 << std::endl;
							std::cout << std::endl;
							
							std::cout << "T =  " << v3 << std::endl;
							std::cout << "V =  " << v4 << std::endl;
							std::cout << "V_c =  " << v5 << std::endl;
							std::cout << std::endl;
							
							myfile_1 << "Autovettore eccitato " << indice_stato_eccitato_c << std::endl;
							myfile_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_1 << "T =  " << v3 << std::endl;
							myfile_1 << "V =  " << v4 << std::endl;
							myfile_1 << "V_c =  " << v5 << std::endl;
							myfile_1 << std::endl;
							
							myfile_misurabili_1 << "Coulomb:" << std::endl;
							myfile_misurabili_1 << "r_i^2 =  " << v1 << std::endl;
							myfile_misurabili_1 << "r_ij^2 =  " << v2 << std::endl;
							myfile_misurabili_1 << "T =  " << v3 << std::endl;
							myfile_misurabili_1 << "V =  " << v4 << std::endl;
							myfile_misurabili_1 << "V_c =  " << v5 << std::endl;
							myfile_misurabili_1 << "V_t =  " << (v4+v5) << std::endl;
							myfile_misurabili_1 << std::endl;
							if(Corpi.M_u() && A==3)	{
								tipo_matrici v6 = valor_medio(EN.MATRICI_NABLA(),autostato);
								
								myfile_misurabili_1 << "E_ppn = " << autovalori_c(indice_stato_eccitato_c) + (m_A + m_B)*v3 + (4./3.)*m_B*v6<< std::endl;
								myfile_misurabili_1 << "E_nnp = " << autovalori_c(indice_stato_eccitato_c) + (m_A - m_B)*v3 - (4./3.)*m_B*v6<< std::endl;
							}
							
							mkl_free_buffers();
						}
						
						std::cerr << "Misurabili  T = " << t5.elapsed() << std::endl;
						myfile << "Misurabili  T = " << t5.elapsed() << std::endl;
					}
					
				}
				
				std::cerr << "K_MAX = " << K_max_temp << " T = " << t.elapsed() << std::endl;
				myfile << "K_MAX = " << K_max_temp << " T = " << t.elapsed() << std::endl;
				
				mkl_free_buffers();
			}
		}
	}
	std::cerr << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile.close();
	myfile_1.close();
	myfile_tabella.close();
	myfile_livelli.close();
	myfile_misurabili.close();
	myfile_tabella_1.close();
	myfile_misurabili_1.close();
	
}

int main (int argc, char *argv[]) {
	
	if (argc == 1){
		std::cerr << "Sintassi: A m_max K_min K_max nucleoni P_clus masse_uguali Parita_Orbitale [modo_H]" << std::endl;
		return true;
	}

	MKLVersion ver;
	
	// Print information on CPU optimization in effect
	
	MKLGetVersion(&ver);
	
	printf("Processor optimization: %s\n",ver.Processor);
	omp_set_num_threads(8);
	mkl_set_dynamic(true);
	
	unsigned short A = ((unsigned short) std::strtoul(argv[1],(char **)0, 10));
	unsigned short m_max = ((unsigned short) std::strtoul(argv[2],(char **)0, 10));
	unsigned short K_min = ((unsigned short) std::strtoul(argv[3],(char **)0, 10));
	unsigned short K_max = ((unsigned short) std::strtoul(argv[4],(char **)0, 10));
	std::vector <char> nucleoni;
	for (short i = 5; i<5+A; ++i) {
		nucleoni.push_back(*argv[i]);
	}
	unsigned short P = ((unsigned short) std::strtoul(argv[5+A],(char **)0, 10));	//0 = pari ; 1 = dispari ; 2 = tutte
	bool controllo = ((bool) std::strtoul(argv[6+A],(char **)0, 10)); //0 = masse diverse ; 1 = masse uguali
	unsigned short Parita_Orbitale = ((unsigned short) std::strtoul(argv[7+A],(char **)0, 10));; //0 = pari ; 1 = dispari ; 2 = tutte
	bool modo_cluster = false;
	beta_riferimento = 1.;
	if (A==4) {
		beta_riferimento = 2.;
		if (argc == 9+A) {
			modo_cluster = ((bool) std::strtoul(argv[8+A],(char **)0, 10)); //mettere a 1 per avere modo H sui 4 corpi
		}
	}
	Init_Corpi_L Corp(nucleoni,controllo);
	std::string Nome_dir_vec = "lista_vettori"; 
	int crea_dir = mkdir(Nome_dir_vec.c_str(), 0777);
	init_main_0(A,m_max,K_min,K_max,Corp,P,modo_cluster,Parita_Orbitale);
	return true;
}


