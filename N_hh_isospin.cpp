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
				 const unsigned short & S_tot,
				 const unsigned short & T_max,
				 const Init_Corpi_LST & Corpi,
				 const bool & modo_clust){
	
	unsigned short ordine = 200;
	unsigned short N = A - 1;
	unsigned short P = 2;
	std::string s, s_0, s_0_p, s_0_d, s_1 , s_l, s_misurabili;
	std::string s_ext_txt = ".txt";
	std::string s_ext_dat = ".dat";
	std::string s_ext_p = "_p";
	std::string s_ext_d = "_d";
	std::string s_cluster,s_masse;
	
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
	
	
	{
		s = "N_HH_" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;
		
		s_0 = "N_HH_tab" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;
		
		s_0_p = "N_HH_tab" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;
		
		s_0_d = "N_HH_tab" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;
		
		s_1 = "N_HH_" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;

		s_l = "N_HH_Livelli_E_T_0_" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;

		s_misurabili = "N_HH_Misurabili_T_" + scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_K_max_[" + scrivi(K_min) + "-" + scrivi(K_max) + "]_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max) ;
	}
	{
		s += s_cluster + s_masse + s_ext_dat;
		
		s_0 += s_cluster + s_masse + s_ext_txt;
		
		s_0_p += s_cluster + s_masse + s_ext_p + s_ext_txt;
		
		s_0_d += s_cluster + s_masse + s_ext_d + s_ext_txt;
		
		s_1 += s_cluster + s_masse + s_ext_txt;
		
		s_l += s_cluster + s_masse + s_ext_txt;
		
		s_misurabili += s_cluster + s_masse + s_ext_txt;
	}
	
	std::string s_vec,s_vec_temp,s_vec_temp_c;
	{
		s_vec = "./lista_vettori/Vec_"+ scrivi(A) + "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max) + "_S_" + scrivi(S_tot) + "_T_" + scrivi(T_max);
		s_vec += s_cluster /*+ s_masse*/;
	}
	std::fstream myfile_vec;
	myfile_vec.precision(10); 		
	bool scrivi_vettore;
	bool scrivi_vettore_c;
	
	std::ofstream myfile;
	std::ofstream myfile_1;
	std::ofstream myfile_tabella;
	std::ofstream myfile_tabella_p;
	std::ofstream myfile_tabella_d;
	std::ofstream myfile_livelli;
	std::ofstream myfile_misurabili;
	{
		myfile.open (s.c_str());
		myfile_1.open (s_1.c_str());
		myfile_tabella.open (s_0.c_str());
		myfile_tabella_p.open (s_0_p.c_str());
		myfile_tabella_d.open (s_0_d.c_str());
		myfile_livelli.open (s_l.c_str());
		myfile_misurabili.open (s_misurabili.c_str());
		myfile.precision(std::numeric_limits<long double>::digits10); 		
		myfile_1.precision(std::numeric_limits<long double>::digits10); 		
		myfile_tabella.precision(10); 		
		myfile_tabella_p.precision(10); 		
		myfile_tabella_d.precision(10);
		myfile_livelli.precision(10);
		myfile_misurabili.precision(10);
		std::cerr.precision(10);
		std::cout.precision(10);
	}
	
	{		
		myfile_1 << "Binding Energy per " << A << " nucleoni : ";
		myfile_tabella << "Binding Energy per " << A << " nucleoni : ";
		myfile_tabella_p << "Binding Energy per " << A << " nucleoni : ";
		myfile_tabella_d << "Binding Energy per " << A << " nucleoni : ";
		myfile_livelli << "Primi Livelli Energetici per " << A << " nucleoni : ";
		for (short i = 0; i<A; ++i) {
			myfile_1<< Corpi.CORPI(i)<< " ";
			myfile_tabella<< Corpi.CORPI(i)<< " ";
			myfile_tabella_p<< Corpi.CORPI(i)<< " ";
			myfile_tabella_d<< Corpi.CORPI(i)<< " ";
			myfile_livelli<< Corpi.CORPI(i)<< " ";
		}
		myfile_1 << std::endl;
		myfile_tabella << std::endl;
		myfile_tabella_p << std::endl;
		myfile_tabella_d << std::endl;
		myfile_livelli << std::endl;
		myfile_1 << " m_max = " << m_max << " K_max = " << K_max << "	S_tot = " << S_tot <<"/2	T_max = " << T_max << "/2" << std::endl;
		myfile_tabella << " m_max = " << m_max << " K_max = " << K_max << "	S_tot = " << S_tot <<"/2	T_max = " << T_max << "/2"<< std::endl;
		myfile_tabella_p << " m_max = " << m_max << " K_max = " << K_max << "	S_tot = " << S_tot <<"/2	T_max = " << T_max << "/2"<< std::endl;
		myfile_tabella_d << " m_max = " << m_max << " K_max = " << K_max << "	S_tot = " << S_tot <<"/2	T_max = " << T_max << "/2"<< std::endl;
		myfile_livelli << " m_max = " << m_max << " K_max = " << K_max << "	S_tot = " << S_tot <<"/2	T_max = " << T_max << "/2"<< std::endl;
	}
	if(Corpi.M_u() && A==3)	{
		myfile_tabella << "K_max	N_hh	Parita'		E		E_m		E_c" << std::endl;
		myfile_tabella_p << "K_max	N_hh	Parita'		E		E_m		E_c" << std::endl;
		myfile_tabella_d << "K_max	N_hh	Parita'		E		E_m		E_c" << std::endl;
	}
	else{
		myfile_tabella << "K_max	N_hh	Parita'		E		E_c" << std::endl;
		myfile_tabella_p << "K_max	N_hh	Parita'		E		E_c" << std::endl;
		myfile_tabella_d << "K_max	N_hh	Parita'		E		E_c" << std::endl;
	}
	
	boost::timer t6;
	unsigned short incremento = 2;
	Gauss_Laguerre G(ordine);
	HyperInt HI(ordine/2);

	Integrazioni_3_P INT_3P(N);
	Accoppiamenti_Spinoriali_S SS (Corpi.D_s(),clus,S_tot,Corpi.D_Z_min());
	Accoppiamenti_Spinoriali_S_max TT (Corpi.D_t(),clus,T_max,Corpi.D_T());
	Accoppiamenti_Iperangolari_K_max KK (N, 0, 0, K_min,P,clus,0);
	
	//	Energia_radiale_isospin ER ( SS, TT,N, G, HI, K_min, m_max, Corpi ,clus);
	//	Energia_isospin EN ( HI, &INT_3P , ER, KK );
	
	//lancsoz
	VETTORI mygen(m_max + 1);
	VETTORI mygen_c(m_max + 1);
	VETTORI mygen1(m_max + 1);
	VETTORI mygen1_c(m_max + 1);
	VETTORI mygen2(m_max + 1);
	VETTORI mygen2_c(m_max + 1);
	typedef boost::lagged_fibonacci607 Gen_i;
	typedef ietl::vectorspace<VETTORI> Vecspace;
	Gen_i mygen_i;
	bool primopasso = true;
	bool primopasso1 = true;
	bool primopasso2 = true;
	bool primopasso_c = true;
	bool primopasso1_c = true;
	bool primopasso2_c = true;	
	tipo_matrici rel_tol = 500000000*std::numeric_limits<double>::epsilon();
	tipo_matrici abs_tol = 1000000*std::pow(std::numeric_limits<double>::epsilon(),2./3);
	int n_lowest_eigenval = 1;	
	//lancsoz

	
	for (unsigned short K_max_temp = K_min; K_max_temp <= K_max; K_max_temp+=incremento) {
		
		//if (K_max_temp == 30) incremento = 10;
		
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
			for (unsigned long i = 0; i < SS.D() ; ++i){
				myfile << i+1 << "	S = [ ";
				for(unsigned short j = 0; j < SS.D_S(i).size() ; ++j){
					myfile << SS.D_S(i)[j];
					myfile << "	";
				}
				myfile << "]	S_z = [ ";
				myfile << SS.D_Z(i);
				myfile << "	";
				myfile << "]";
				myfile << "	parita = [ ";
				for(unsigned short j = 0; j < SS.P_S(i).size() ; ++j){
					myfile << SS.P_S(i)[j];
					myfile << "	";
				}
				myfile << "]	 somma quadra dei coefficienti dello sviluppo = ";
				myfile << quad_sum(SS.cm(i));				
				myfile << std::endl;
			}
			myfile << "Dimensione della base - SPIN = " << SS.D() << " combinazioni." << std::endl;
			for (unsigned long i = 0; i < TT.D() ; ++i){
				myfile << i+1 << "	T = [ ";
				for(unsigned short j = 0; j < TT.D_S(i).size() ; ++j){
					myfile << TT.D_S(i)[j];
					myfile << "	";
				}
				myfile << "]	T_z = [ ";
				myfile << TT.D_Z(i);
				myfile << "	";
				myfile << "]";
				myfile << "	parita = [ ";
				for(unsigned short j = 0; j < TT.P_S(i).size() ; ++j){
					myfile << TT.P_S(i)[j];
					myfile << "	";
				}
				myfile << "]	 somma quadra dei coefficienti dello sviluppo = ";
				myfile << quad_sum(TT.cm(i));				
				myfile << std::endl;
			}
			myfile << "Dimensione della base - ISOSPIN = " << TT.D() << " combinazioni." << std::endl;
			
		}
		boost::timer t1;
		{
			
			std::cerr << "Energia  T = " << t1.elapsed() << std::endl;
			myfile << "Energia  T = " << t1.elapsed() << std::endl;
			
		}
		Accoppiamenti_Parita BASE(KK,SS,TT,m_max + 1);
		int max_iter = 5*BASE.D();

		std::cerr << "BASE.D() "<< BASE.D() << std::endl;
		myfile_tabella_p << K_max_temp << "	";
		myfile_tabella_d << K_max_temp << "	";
		myfile_tabella_p << BASE.D_h_p() << std::endl;
		myfile_tabella_d << BASE.D_h_d() << std::endl;
	}
	std::cerr << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile << "Esecuzione  T = " << t6.elapsed() << std::endl;
	myfile.close();
	myfile_1.close();
	myfile_tabella.close();
	myfile_tabella_p.close();
	myfile_tabella_d.close();
	myfile_livelli.close();
	myfile_misurabili.close();
}

int main (int argc, char *argv[]) {
	
	if (argc == 1){
		std::cerr << "Sintassi: A m_max K_min K_max nucleoni S_tot T_max masse_uguali [modo_H]" << std::endl;
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
	unsigned short S_tot = ((unsigned short) std::strtoul(argv[5+A],(char **)0, 10));
	unsigned short T_max = ((unsigned short) std::strtoul(argv[6+A],(char **)0, 10));
	bool controllo = ((bool) std::strtoul(argv[7+A],(char **)0, 10));
	bool modo_cluster = false;
	if (A==4) {
		beta_riferimento = 2.;
		if (argc == 9+A) {
			modo_cluster = ((bool) std::strtoul(argv[8+A],(char **)0, 10)); //mettere a 1 per avere modo H sui 4 corpi
		}				
	}
	
	Init_Corpi_LST Corp(nucleoni,controllo);
	std::string Nome_dir_vec = "lista_vettori"; 
	int crea_dir = mkdir(Nome_dir_vec.c_str(), 0777);
	init_main_0(A,m_max,K_min,K_max,S_tot,T_max,Corp,modo_cluster);
	return true;
}


