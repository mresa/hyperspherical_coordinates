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
				 const unsigned short & K_max,
				 const Init_Corpi_L & Corpi,
				 const unsigned short & P,
				 const bool & modo_clust,
				 const unsigned short & P_orb){
	
	unsigned short ordine = 200;
	unsigned short N = A - 1;
	
	std::string s_f_y, s_f;
	std::string s_k,s_masse;
	std::string s_ext_txt = "_ecc.txt";
	std::string s_corpi = "_" + Corpi.scrivi_corpi() + "__m_max_" + scrivi(m_max);
	std::string s_deut = "_deut";
	
	std::vector <unsigned short> clus;
	
	clus.push_back(A);
	if(Corpi.M_u()){
		s_masse = "_M_u";
	}
	else {
		s_masse = "";
	}
	
	Accoppiamenti_Iperangolari_K_max KK (N, 0, 0, K_max,P,clus,P_orb);
	s_k = KK.filename();
	{
		s_f_y = "f_quad_";
		s_f = "f_tot_";
	}
	{
		s_f_y += s_corpi + s_k + s_masse + s_ext_txt;
		s_f += s_corpi + s_k + s_masse + s_ext_txt;
	}
	std::string s_deut,s_stato,s_stato_c;
	{
		s_deut = "./lista_vettori/";
		s_stato = "./lista_vettori/";
		s_stato_c = "./lista_vettori/";
	}
	std::fstream myfile_vec;
	myfile_vec.precision(10); 		
	
	std::ofstream myfile_f_y;
	std::ofstream myfile_f;
	{
		myfile_f_y.open (s_f_y.c_str());
		myfile_f.open (s_f.c_str());
		myfile_f_y.precision(std::numeric_limits<long double>::digits10); 		
		myfile_f.precision(std::numeric_limits<long double>::digits10); 		
		std::cerr.precision(10);
		std::cout.precision(10);
	}
	
	Gauss_Laguerre G(ordine);
	VETTORI v_deutone(51);
	VETTORI v_stato((m_max+1)*KK.D());
	VETTORI v_stato_c((m_max+1)*KK.D());
	
	//carica deutone
	myfile_vec.open(s_deut.c_str(),std::ios::in);
	if (myfile_vec) {
		myfile_vec >> v_deutone;
	}
	myfile_vec.close();

	//carica stato
	myfile_vec.open(s_stato.c_str(),std::ios::in);
	if (myfile_vec) {
		myfile_vec >> v_stato;
	}
	myfile_vec.close();

	//carica stato_c
	myfile_vec.open(s_stato_c.c_str(),std::ios::in);
	if (myfile_vec) {
		myfile_vec >> v_stato_c;
	}
	myfile_vec.close();

	//costruisci la matrice A_mn
	MATRICI Coeff_mn(A_mn(KK,m_max+1,v_stato));
	MATRICI Coeff_mn_c(A_mn(KK,m_max+1,v_stato));
	
	std::vector<long double> funz_y;
	std::vector<long double> funz_y_c;
	long double funz_y_tot=0.;
	long double funz_y_c_tot=0.;
	for (int ky = 0; ky < G.N() ; ++ky){ //ipsilon  *G.X(ky)*G.X(ky)
		long double funz_y_temp=0.;
		long double funz_y_c_temp=0.;
		for (int k = 0; k < G.N() ; ++k)//x
			for (int m = 0; m < Coeff_mn.size1() ; ++m)
				for (int n = 0; n < Coeff_mn.size2() ; ++n){
					long double temp=0.;
#pragma omp parallel for reduction(+:temp) schedule( guided )
					for (int m_d = 0; m_d < v_deutone.size() ; ++m_d){
						temp += costanti_integrando_deutone(beta,m,n,m_d)*G.W(k)*integrando_deutone(beta,m,n,m_d,G.X(k),G.X(ky));
					}
					funz_y_temp +=Coeff_mn(m,n)*temp;
					funz_y_c_temp +=Coeff_mn_c(m,n)*temp;
				}
		funz_y.push_back(pow(funz_y_temp*G.X(ky)*G.X(ky),2));
		funz_y_c.push_back(pow(funz_y_c_temp*G.X(ky)*G.X(ky),2));
		//scrivi i parziali nei file
		funz_y_tot +=funz_y.back()*exp (G.X(ky))*G.W(ky);
		funz_y_c_tot +=funz_y_c.back()*exp (G.X(ky))*G.W(ky);
	}
	//scrivi i totali nei file

	myfile_f_y.close();
	myfile_f.close();
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


