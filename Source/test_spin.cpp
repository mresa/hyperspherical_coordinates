/*
 *  test_spin.cpp
 *  Tesi
 *
 *  Created by Marco on 15/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

/*
 *  Test_Spin.cpp
 *  Tesi
 *
 *  Created by Marco Resa on 08/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


/*
 #include <algorithm>
 #include <gsl/gsl_sf_coupling.h>
 #include "algebra.h"
 #include "cout_vec.h"
 #include <time.h>
 #include <stdio.h>
 #include "triangolari.h"
 #include "accoppiamenti_spinoriali.h"
 #include "err.h"
 
 
 int main(int argc, char *argv[]){
 std::cout.precision(22);	
 using namespace std;
 
 // unsigned short s1 = ((unsigned short) strtoul(argv[1],(char **)0, 10));
 // unsigned short s2 = ((unsigned short) strtoul(argv[2],(char **)0, 10));
 // unsigned short s3 = ((unsigned short) strtoul(argv[3],(char **)0, 10));
 // unsigned short S1 = ((unsigned short) strtoul(argv[4],(char **)0, 10));
 // unsigned short S2 = ((unsigned short) strtoul(argv[5],(char **)0, 10));
 // short Z = ((unsigned short) strtoul(argv[6],(char **)0, 10));
 
 unsigned short s1 = 1;
 unsigned short s2 = 1;
 unsigned short s3 = 1;
 //unsigned short S1 = ((unsigned short) strtoul(argv[1],(char **)0, 10));
 std::vector <unsigned short> s;
 s.push_back(s1);
 s.push_back(s2);
 s.push_back(s3);
 
 //std::vector <unsigned short> S;
 //S.push_back(S1);
 //	S.push_back(S2);
 
 Accoppiamenti_Spinoriali_S ACC(s);
 
 cout << endl;
 cout << endl;
 cout << endl;
 std::cout << " **************				POST PROGRAMMA				**************"<< std::endl;
 
 cout << endl;
 cout << endl;
 cout << endl;
 cout << ACC.D() << " stati  ";
 cout << endl;
 
 
 for(unsigned short stato = 0 ; stato< ACC.D(); ++stato){
 cout << " s = ";
 cout_vec(ACC.D_s());
 cout << " S = ";
 cout_vec(ACC.D_S(stato));
 cout << " Z = " << ACC.D_Z(stato);
 cout << endl;
 cout << " Sviluppo in ";
 cout << ACC.Dz(stato);
 cout << " stati ";
 cout << endl;
 
 for(unsigned short sviluppo = 0 ; sviluppo< ACC.Dz(stato); ++sviluppo){
 cout_vec(ACC.D_z(stato,sviluppo));
 cout << "		" << ACC.cm(stato,sviluppo) << endl;
 }
 cout << "Somma Quadra dei coefficienti = " << quad_sum(ACC.cm(stato)) << endl;
 if ( std::abs(quad_sum(ACC.cm(stato)) - 1) > 1000*err) {
 cout << endl;
 cout << endl;
 cout << "					ERRORE !!! " << endl;
 cout << endl;
 cout << endl;
 
 }
 }
 
 return true;
 
 }
 */

#include <algorithm>
#include <gsl/gsl_sf_coupling.h>
#include "algebra.h"
#include "cout_vec.h"
#include <time.h>
#include <stdio.h>
#include "triangolari.h"
#include "accoppiamenti_spinoriali.h"
#include "err.h"


int main(int argc, char *argv[]){
	std::cout.precision(22);	
	using namespace std;
	
	unsigned short s1 = 1;
	unsigned short s2 = 1;
	unsigned short s3 = 1;
	unsigned short s4 = 1;
	//	unsigned short S = ((unsigned short) strtoul(argv[1],(char **)0, 10));
	bool cl = ((bool) strtoul(argv[1],(char **)0, 10));
	std::vector <unsigned short> s;
	s.push_back(s1);
	s.push_back(s2);
	s.push_back(s3);
	s.push_back(s4);
	std::vector < unsigned short >  CLUSTER;
	if(cl){
		CLUSTER.push_back(2);
		CLUSTER.push_back(2);
	}
	else{
		CLUSTER.push_back(4);
	}
	
	Accoppiamenti_Spinoriali_S_max ACC(s,CLUSTER,2,0);
	
	cout << endl;
	cout << endl;
	cout << endl;
	std::cout << " **************				POST PROGRAMMA				**************"<< std::endl;
	
	cout << endl;
	cout << endl;
	cout << endl;
	cout << ACC.D() << " stati  ";
	cout << endl;
	
	
	for(unsigned short stato = 0 ; stato< ACC.D(); ++stato){
		cout << " s = ";
		cout_vec(ACC.D_s());
		cout << " S = ";
		cout_vec(ACC.D_S(stato));
		cout << " Z = " << ACC.D_Z(stato);
		cout << endl;
		cout << " Parita = ";
		cout_vec(ACC.P_S(stato));
		cout << endl;
		cout << " Sviluppo in ";
		cout << ACC.Dz(stato);
		cout << " stati ";
		cout << endl;
		
		for(unsigned short sviluppo = 0 ; sviluppo< ACC.Dz(stato); ++sviluppo){
			cout_vec(ACC.D_z(stato,sviluppo));
			cout << "		" << ACC.cm(stato,sviluppo) << endl;
		}
		cout << "Somma Quadra dei coefficienti = " << quad_sum(ACC.cm(stato)) << endl;
		if ( std::abs(quad_sum(ACC.cm(stato)) - 1) > 1000*err) {
			cout << endl;
			cout << endl;
			cout << "					ERRORE !!! " << endl;
			cout << endl;
			cout << endl;
			
		}
	}
	
	return true;
	
}