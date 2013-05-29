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
#include "proiezione_deutone.h"


int main (int argc, char *argv[]) {
	MKLVersion ver;
	MKLGetVersion(&ver);
	
	printf("Processor optimization: %s\n",ver.Processor);
	omp_set_num_threads(8);
	mkl_set_dynamic(true);
	
	unsigned short A = 3;
	unsigned short ordine = 200;
	unsigned short N = A - 1;
	unsigned short m_max = 30;
	unsigned short K_max = 60;
	unsigned short P = 0;
	unsigned short P_orb = 0;
	
	std::vector <unsigned short> clus;
	clus.push_back(A);

	Accoppiamenti_Iperangolari_K_max KK (N, 0, 0, K_max,P,clus,P_orb);
	std::cerr << "Accoppiamenti_Iperangolari_K_max creati" << std::endl;

	Gauss_Laguerre G(ordine);

	std::string s_vec_fond,s_vec_ecc,s_mat_fond,s_mat_ecc;
	s_vec_fond = "./deutone/Vec_fond.txt";
	s_vec_ecc = "./deutone/Vec_ecc.txt";
	s_mat_fond = "./deutone/Mat_fond.txt";
	s_mat_ecc = "./deutone/Mat_ecc.txt";

	std::fstream myfile;
	myfile.precision(10); 		

	VETTORI vettore_caricato;

	myfile.open(s_vec_fond.c_str(),std::ios::in);
	myfile >> vettore_caricato;
	myfile.close();
	std::cerr << "vettore caricato : dimensione = " << vettore_caricato.size() << std::endl;
	MATRICI MAT(A_mn(KK,m_max+1,vettore_caricato));
	std::cerr << "matrice creata : dimensione = [" << MAT.size1() << ","<< MAT.size2() << "]" << std::endl;

	myfile.open(s_mat_fond.c_str(),std::ios::out);
	myfile << MAT << std::endl;
	myfile.close();
	std::cerr << "matrice scritta" << std::endl;

	myfile.open(s_vec_ecc.c_str(),std::ios::in);
	myfile >> vettore_caricato;
	myfile.close();
	std::cerr << "vettore caricato : dimensione = " << vettore_caricato.size() << std::endl;

	MAT = A_mn(KK,m_max+1,vettore_caricato);
	std::cerr << "matrice creata : dimensione = [" << MAT.size1() << ","<< MAT.size2() << "]" << std::endl;

	myfile.open(s_mat_ecc.c_str(),std::ios::out);
	myfile << MAT << std::endl;
	myfile.close();
	std::cerr << "matrice scritta" << std::endl;
	return true;
}


