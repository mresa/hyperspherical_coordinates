#include "preambolo.h"
#include "gauss-laguerre.h"
#include <boost/timer.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <assert.h>
#include <boost/random.hpp>
#include <boost/limits.hpp>
#include <cmath>
#include <limits>
#include "proiezione_deutone.h"


int main (int argc, char *argv[]) {
	MKLVersion ver;
	MKLGetVersion(&ver);
	
	printf("Processor optimization: %s\n",ver.Processor);
	omp_set_num_threads(8);
	mkl_set_dynamic(true);
	
	unsigned short ordine = 300;
	Gauss_Laguerre G(ordine);

	std::string s,s_p_y,s_p,s_vec_deutone,s_matrice,s_a_y,s_raggioq;
	s = "./deutone/fond.txt";
	s_raggioq = "./deutone/rq_fond.txt";
	s_p_y = "./deutone/P_y_fond.txt";
	s_a_y = "./deutone/A_y_fond.txt";
	s_p = "./deutone/P_fond.txt";
	s_vec_deutone = "./deutone/Vec_d_fond.txt";
	s_matrice = "./deutone/Mat_fond.txt";

	long double beta = 3.;

	
	/*******************************************************************/
	
	std::fstream myfile;
	myfile.precision(10); 		
	std::fstream myfile_s;
	myfile_s.precision(10); 		
	myfile_s.open(s.c_str(),std::ios::out);
	
	myfile.open(s_raggioq.c_str(),std::ios::out);
	myfile.close();
	myfile.open(s_a_y.c_str(),std::ios::out);
	myfile.close();
	myfile.open(s_p_y.c_str(),std::ios::out);
	myfile.close();
	myfile.open(s_p.c_str(),std::ios::out);
	myfile.close();
	
	VETTORI vettore_caricato;
	MATRICI matrice_caricata;
	
	myfile.open(s_vec_deutone.c_str(),std::ios::in);
	myfile >> vettore_caricato;
	myfile.close();
	std::cerr << "vettore caricato : dimensione = " << vettore_caricato.size() << std::endl;
	std::cerr << "somma quadra vettore  = " << quad_sum_vec(vettore_caricato) << std::endl;
	myfile_s << "vettore caricato : dimensione = " << vettore_caricato.size() << std::endl;
	myfile_s << "somma quadra vettore  = " << quad_sum_vec(vettore_caricato) << std::endl;
	
	myfile.open(s_matrice.c_str(),std::ios::in);
	myfile >> matrice_caricata;
	myfile.close();
	
	std::cerr << "matrice caricata : dimensione = [" << matrice_caricata.size1() << ","<< matrice_caricata.size2() << "]" << std::endl;
	std::cerr << "somma quadra matrice  = " << quad_sum_mat(matrice_caricata) << std::endl;
	
	myfile_s << "matrice caricata : dimensione = [" << matrice_caricata.size1() << ","<< matrice_caricata.size2() << "]" << std::endl;
	myfile_s << "somma quadra matrice  = " << quad_sum_mat(matrice_caricata) << std::endl;
	
	long double probabilita_totale = 0.;
	long double raggio_q = 0.;
	long double probabilita_y,incremento,incremento_rq;
	std::cerr << "k = " ;
	myfile_s << "k = " ;
	for (int k = 0; k < G.N() ; ++k){
		std::cerr <<  k << " "<< std::endl;
		myfile_s <<  k << " "<< std::endl;
		
		probabilita_y = integrale_deutone_y(G,beta,vettore_caricato,matrice_caricata,G.X(k));
		
		myfile.open(s_a_y.c_str(),std::ios::out | std::ios::app);
		myfile << G.X(k) << " " << probabilita_y << std::endl;
		myfile.close();
		std::cerr << "Ampiezza: " << G.X(k) << " " << probabilita_y << std::endl;
		myfile_s << "Ampiezza: " << G.X(k) << " " << probabilita_y << std::endl;
		
		probabilita_y *= probabilita_y;
		
		myfile.open(s_p_y.c_str(),std::ios::out | std::ios::app);
		myfile << G.X(k) << " " << probabilita_y << std::endl;
		myfile.close();
		std::cerr << "probabilita_y: " << G.X(k) << " " << probabilita_y << std::endl;
		myfile_s << "probabilita_y: " << G.X(k) << " " << probabilita_y << std::endl;
		
		incremento_rq = probabilita_y	*	pow (G.X(k),4)	*	G.W(k)	*	exp(G.X(k))	;
		incremento = probabilita_y	*	pow (G.X(k),2)	*	G.W(k)	*	exp(G.X(k))	;
		
		std::cerr << "incremento: " << G.X(k) << " " << incremento << std::endl;
		myfile_s << "incremento: " << G.X(k) << " " << incremento << std::endl;
		std::cerr << "incremento rq: " << G.X(k) << " " << incremento << std::endl;
		myfile_s << "incremento rq: " << G.X(k) << " " << incremento << std::endl;
		
		if (incremento != incremento ) {
			std::cerr << "errore in main" << std::endl;
			myfile_s << "errore in main" << std::endl;
			incremento = 0.;
		}
		if (incremento_rq != incremento_rq ) {
			std::cerr << "errore in main rq" << std::endl;
			myfile_s << "errore in main rq" << std::endl;
			incremento_rq = 0.;
		}
		probabilita_totale += incremento;
		raggio_q += incremento_rq;
		std::cerr << "probabilita_totale = " << probabilita_totale << std::endl;
		myfile_s << "probabilita_totale = " << probabilita_totale << std::endl;
		myfile.open(s_p.c_str(),std::ios::out | std::ios::app);
		myfile << "probabilita_totale_parziale = " << probabilita_totale << std::endl;
		myfile.close();					
		myfile.open(s_raggioq.c_str(),std::ios::out | std::ios::app);
		myfile << "rq parziale = " << raggio_q << std::endl;
		myfile.close();					
	}
	std::cerr << std::endl;
	std::cerr << "probabilita_totale = " << probabilita_totale << std::endl;
	myfile_s << std::endl;
	myfile_s << "probabilita_totale = " << probabilita_totale << std::endl;
	
	myfile.open(s_p.c_str(),std::ios::out | std::ios::app);
	myfile << probabilita_totale << std::endl;
	myfile.close();
	
	myfile.open(s_raggioq.c_str(),std::ios::out | std::ios::app);
	myfile << raggio_q/probabilita_totale << std::endl;
	myfile.close();					
	
	std::cerr << "probabilita_totale scritta " << std::endl;
	myfile_s << "probabilita_totale scritta " << std::endl;
	myfile_s.close();
	return true;
}

