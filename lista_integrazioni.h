/*
 *  lista_integrazioni.h
 *  Tesi
 *
 *  Created by Marco on 30/04/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef LISTA_INTEGRAZIONI_H
#define LISTA_INTEGRAZIONI_H

#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include "tipi_definiti.h"
#include "hyper.h"
#include "hyperintegrazione.h"
#include "io_sparse.h"
#include <iostream>
#include <ext/hash_map>
#include <string.h>
#include "err.h"

inline long double hyper_int_ripelle (const HyperInt & G,
									  const unsigned short i,
									  std::vector <unsigned short> P_1,
									  std::vector <unsigned short> P_2,
									  std::vector <unsigned short> P_3){
	unsigned short ordine = G.N();
	
	assert (ordine != 0);
	
	if (i == 0) return 1.;
	
	else if (dis_triang(P_2[0], P_1[0], P_3[0]) == false || dis_triang(P_2[2], P_1[2], P_3[2]) == false ) return 0.;
	
	
	else { 
		long double integrale = 0.;
		
		unsigned short k1[3] =  { (P_1[0] - 2*P_1[1] - P_1[2]) , (P_2[0] - 2*P_2[1] - P_2[2]) , (P_3[0] - 2*P_3[1] - P_3[2]) } ;
		unsigned short dub_alpha = k1[0] + k1[1] + k1[2] + 3*i - 2; //attenzione: sarebbe 3i -5 ma il mio i parte da 0
		unsigned short dub_beta = P_1[2] + P_2[2] + P_3[2] + 1;
		unsigned short n[3] =  { P_1[1] , P_2[1] , P_3[1] };
		long double alpha[3] =  { k1[0] + 3. * i / 2. - 1. , k1[1] + 3. * i /2. - 1. , k1[2] + 3. * i /2. - 1. };
		long double beta[3] =  { P_1[2] + 1./2. , P_2[2] + 1./2. , P_3[2] + 1./2.};
		
		if(dub_alpha %2 == 0){ //pari pari -> gauss
			if(dub_beta %2 == 0){
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short I = 0 ; I < ordine; ++I) {
					long double temp = 1.;				
					for (short j = 0 ; j < 3; ++j) {
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  sqrt(pow(1. - G.X(i,dub_alpha,dub_beta) , dub_alpha))	* 
								  sqrt(pow(1. + G.X(i,dub_alpha,dub_beta) , dub_beta))	* 
								  G.W(i,dub_alpha,dub_beta));
				}
			}
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					
					long double temp = 1.;
					
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= (sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *
								 jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta)));
					}
					
					integrale += (temp *
								  sqrt(pow(1. - G.X(i,dub_alpha,dub_beta) , dub_alpha))	* 
								  sqrt(pow(1. + G.X(i,dub_alpha,dub_beta) , dub_beta-1))	* 
								  G.W(i,dub_alpha,dub_beta));
				}				
			}
		}
		
		else {
			// alpha semintero, beta intero
			if(dub_beta %2 == 0){  // legendre 2 non modificata
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2 )) *
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
					}
					
					integrale += (temp *
								  sqrt(pow(1. - G.X(i,dub_alpha,dub_beta) , dub_alpha-1))	* 
								  sqrt(pow(1. + G.X(i,dub_alpha,dub_beta) , dub_beta))	* 
								  G.W(i,dub_alpha,dub_beta));
				}
			}
			// alpha semintero, beta semintero
			else{
#pragma omp parallel for reduction(+:integrale) schedule( guided )
				for (short i = 0 ; i < ordine; ++i) {
					long double temp = 1.;
					for (short j = 0 ; j < 3; ++j) {
						
						temp *= sqrt(pow ((long double) 2., alpha[j] + beta[j] + 2)) *						
						jacobi_N( n[j], alpha[j], beta[j], G.X(i,dub_alpha,dub_beta));
						
					}
					integrale += (temp *
								  sqrt(pow(1. - G.X(i,dub_alpha,dub_beta) , dub_alpha-1))	*
								  sqrt(pow(1. + G.X(i,dub_alpha,dub_beta) , dub_beta-1))	*
								  G.W(i,dub_alpha,dub_beta));
				}
			}
		}
		integrale /= sqrt(pow ((long double) 2.,dub_alpha + dub_beta + 4));
		return integrale;

	}
}

class Numeri_3_P { //simmetrico per scambio di tripletti, con stesso i
	unsigned short indice_N;
	std::vector <unsigned short> P_1;
	std::vector <unsigned short> P_2;
	std::vector <unsigned short> P_3;
	unsigned long ValHash;
	MATRICI_INTERE M;
public:
	Numeri_3_P(const unsigned short & Indice,
			   const MATRICI_INTERE & M_i) {
		indice_N = Indice;
		if (indice_N == 0) {
			ValHash=0.;
		}
		else {
			M=M_i;
			for(short i = 0; i<3;++i){
				P_1.push_back(M(i,0));
				P_2.push_back(M(i,1));
				P_3.push_back(M(i,2));
			}
			
			//deve essere simmetrico per scambi 1-2-3 ma non per altri scambi
			ValHash = P_1[0]*P_1[1]*P_1[2] + P_2[0]*P_2[1]*P_2[2] + P_3[0]*P_3[1]*P_3[2]; 
			ValHash += 3*( P_1[0]*P_2[0]*P_3[0] + P_1[1]*P_2[1]*P_3[1] + P_1[2]*P_2[2]*P_3[2]); 
			ValHash += 5*( (P_1[0] + P_2[0] +P_3[0]) * (P_1[0] + P_2[0] +P_3[0]) );
			ValHash += 5*( (P_1[1] + P_2[1] +P_3[1]) * (P_1[1] + P_2[1] +P_3[1]) );
			ValHash += 5*( (P_1[2] + P_2[2] +P_3[2]) * (P_1[2] + P_2[2] +P_3[2]) );
			ValHash += 7*( (P_1[0] + P_1[1] +P_1[2]) * (P_1[0] + P_1[1] +P_1[2]) );
			ValHash += 7*( (P_2[0] + P_2[1] +P_2[2]) * (P_2[0] + P_2[1] +P_2[2]) );
			ValHash += 7*( (P_3[0] + P_3[1] +P_3[2]) * (P_3[0] + P_3[1] +P_3[2]) );
			ValHash += 11* (P_1[0] + P_1[1] + P_1[2] + P_2[0] + P_2[1] + P_2[2] + P_3[0] + P_3[1] + P_3[2]);			
		}
	}
	Numeri_3_P(const unsigned short & Indice,
			   const hyper & H1,
			   const hyper & H2,
			   const hyper & H3) {
		indice_N = Indice;
		if (indice_N == 0) {
			ValHash=0.;
		}
		else {
			{
				P_1.push_back(H1.k(indice_N));
				P_1.push_back(H1.n(indice_N));
				P_1.push_back(H1.l(indice_N));
			}
			{
				P_2.push_back(H2.k(indice_N));
				P_2.push_back(H2.n(indice_N));
				P_2.push_back(H2.l(indice_N));
			}
			{
				P_3.push_back(H3.k(indice_N));
				P_3.push_back(H3.n(indice_N));
				P_3.push_back(H3.l(indice_N));
			}
			MATRICI_SPARSE M_t(3,3);
			for(short i = 0; i<3;++i){
				M_t(i,0)=P_1[i];
				M_t(i,1)=P_2[i];
				M_t(i,2)=P_3[i];
			}
			M=M_t;
			//deve essere simmetrico per scambi 1-2-3 ma non per altri scambi
			ValHash = P_1[0]*P_1[1]*P_1[2] + P_2[0]*P_2[1]*P_2[2] + P_3[0]*P_3[1]*P_3[2]; 
			ValHash += 3*( P_1[0]*P_2[0]*P_3[0] + P_1[1]*P_2[1]*P_3[1] + P_1[2]*P_2[2]*P_3[2]); 
			ValHash += 5*( (P_1[0] + P_2[0] +P_3[0]) * (P_1[0] + P_2[0] +P_3[0]) );
			ValHash += 5*( (P_1[1] + P_2[1] +P_3[1]) * (P_1[1] + P_2[1] +P_3[1]) );
			ValHash += 5*( (P_1[2] + P_2[2] +P_3[2]) * (P_1[2] + P_2[2] +P_3[2]) );
			ValHash += 7*( (P_1[0] + P_1[1] +P_1[2]) * (P_1[0] + P_1[1] +P_1[2]) );
			ValHash += 7*( (P_2[0] + P_2[1] +P_2[2]) * (P_2[0] + P_2[1] +P_2[2]) );
			ValHash += 7*( (P_3[0] + P_3[1] +P_3[2]) * (P_3[0] + P_3[1] +P_3[2]) );
			ValHash += 11* (P_1[0] + P_1[1] + P_1[2] + P_2[0] + P_2[1] + P_2[2] + P_3[0] + P_3[1] + P_3[2]);			
		}
	}
		
	inline unsigned short I() const{return indice_N;}
	inline unsigned long ValoreHash() const{return ValHash;}
	inline MATRICI_INTERE Mat() const { return M;}
	inline std::vector <unsigned short> P1() const {return P_1;}
	inline std::vector <unsigned short> P2() const {return P_2;}
	inline std::vector <unsigned short> P3() const {return P_3;}

	
};
class Hash_Numeri_3_P {
public:
	unsigned long operator()(const Numeri_3_P & N1) const {
		return N1.ValoreHash();
	}
};
class Compara_Numeri_3_P {
public:
	bool operator()(const Numeri_3_P & N1 ,const Numeri_3_P & N2) {
		if (N1.I()!= N2.I()) return false;
		else if (delta_vec(N1.P1(),N2.P1())){
			if (delta_vec(N1.P2(),N2.P2())){
				if (delta_vec(N1.P3(),N2.P3())) return true;
				else return false;
				
			}
			else if (delta_vec(N1.P2(),N2.P3())){
				if (delta_vec(N1.P3(),N2.P2())) return true;
				else return false;
			}
			else return false;
			
		}
		else if (delta_vec(N1.P1(),N2.P2())){
			if (delta_vec(N1.P2(),N2.P1())){
				if (delta_vec(N1.P3(),N2.P3())) return true;
				else return false;
			}
			else if (delta_vec(N1.P2(),N2.P3())){
				if (delta_vec(N1.P3(),N2.P1())) return true;
				else return false;
			}
			else return false;
			
		}
		else if (delta_vec(N1.P1(),N2.P3())){
			if (delta_vec(N1.P2(),N2.P2())){
				if (delta_vec(N1.P3(),N2.P1())) return true;
				else return false;
			}
			else if (delta_vec(N1.P2(),N2.P1())){
				if (delta_vec(N1.P3(),N2.P2())) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	} 
};
long double Valori_3_P (const HyperInt & G, const Numeri_3_P& N1){return hyper_int_ripelle(G,N1.I(),N1.P1(),N1.P2(),N1.P3());}
typedef __gnu_cxx::hash_map < Numeri_3_P, long double , Hash_Numeri_3_P, Compara_Numeri_3_P > Mappe_3_P;


class Numeri_3_L { //simmetrico per scambio 2/3
	unsigned short L_1;
	unsigned short L_2;
	unsigned short L_3;
	unsigned long ValHash;
	MATRICI_INTERE M;
public:
	Numeri_3_L(const MATRICI_INTERE & M_i) { 
		M=M_i;
		L_1=M(0,0);
		L_2=M(0,1);
		L_3=M(0,2);
		
		ValHash = 11*(L_1 + L_2 + L_3);
		ValHash += 7*(L_3*(L_1 + L_2));
		ValHash += 5*(L_1 + L_2 + L_3)*(L_1 + L_2 + L_3);
		ValHash += L_1*L_2*L_3;
	}
	Numeri_3_L(const unsigned short & Indice,
			   const hyper & H1,
			   const hyper & H2,
			   const hyper & H3) { //ordinate <H1,H2|H3>
		L_1 = H1.l(Indice);
		L_2 = H2.l(Indice);
		L_3 = H3.l(Indice);
		ValHash = 11*(L_1 + L_2 + L_3);
		ValHash += 7*(L_3*(L_1 + L_2));
		ValHash += 5*(L_1 + L_2 + L_3)*(L_1 + L_2 + L_3);
		ValHash += L_1*L_2*L_3;
		
		MATRICI_INTERE M_t(1,3);
		M_t(0,0)=L_1;
		M_t(0,1)=L_2;
		M_t(0,2)=L_3;
		M=M_t;
		
	}
	
	inline unsigned long ValoreHash() const{return ValHash;}
	inline MATRICI_INTERE Mat() const {return M;}
	inline unsigned short L1() const {return L_1;}
	inline unsigned short L2() const {return L_2;}
	inline unsigned short L3() const {return L_3;}
	
};
class Hash_Numeri_3_L {
public:
	unsigned long operator()(const Numeri_3_L & N1) const {
		return N1.ValoreHash();
	}
};
class Compara_Numeri_3_L {
public:
	bool operator()(const Numeri_3_L & N1 ,const Numeri_3_L & N2) {
		if (N1.L3() == N2.L3()){
			if (N1.L1() == N2.L1()){
				if (N1.L2() == N2.L2()) return true;
				else return false;
			}
			else if (N1.L1() == N2.L2()){
				if (N1.L2() == N2.L1()) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	} 
};
long double Valori_3_L (const Numeri_3_L& N1){//ordinate <H1,H2|H3>
	return int_3j_ridotto (N1.L1(),  N1.L2() , N1.L3());
}
typedef __gnu_cxx::hash_map<Numeri_3_L, long double , Hash_Numeri_3_L, Compara_Numeri_3_L> Mappe_3_L;


class Numeri_3_M { //simmetrico per scambio 2/3
	unsigned short L_1;
	unsigned short L_2;
	unsigned short L_3;
	short M_1;
	short M_2;
	short M_3;
	unsigned long ValHash;
	MATRICI_INTERE M;
public:
	Numeri_3_M(const MATRICI_INTERE & M_i) { //ordinate <H1,H2|H3>
		M=M_i;
		L_1 = M(0,0);
		L_2 = M(0,1);
		L_3 = M(0,2);
		M_1 = M(1,0);
		M_2 = M(1,1);
		M_3 = M(1,2);
		
		ValHash = 19*(L_1 + L_2 + L_3);
		ValHash += 7*(L_3*(L_1 + L_2));
		ValHash += 5*(L_1 + L_2 + L_3)*(L_1 + L_2 + L_3);
		ValHash += L_1*L_2*L_3;
		
		ValHash += 19*std::abs(M_1 + M_2 + M_3);
		ValHash += 17*std::abs(M_3*(M_1 + M_2));
		ValHash += 5*(M_1 + M_2 + M_3)*(M_1 + M_2 + M_3);
		ValHash += 3*std::abs(M_1*M_2*M_3);
		
		
	}
	Numeri_3_M(const unsigned short & Indice,
			   const hyper & H1,
			   const hyper & H2,
			   const hyper & H3) { //ordinate <H1,H2|H3>
		L_1 = H1.l(Indice);
		L_2 = H2.l(Indice);
		L_3 = H3.l(Indice);
		M_1 = H1.m(Indice);
		M_2 = H2.m(Indice);
		M_3 = H3.m(Indice);
		
		ValHash = 19*(L_1 + L_2 + L_3);
		ValHash += 7*(L_3*(L_1 + L_2));
		ValHash += 5*(L_1 + L_2 + L_3)*(L_1 + L_2 + L_3);
		ValHash += L_1*L_2*L_3;
		
		ValHash += 19*std::abs(M_1 + M_2 + M_3);
		ValHash += 17*std::abs(M_3*(M_1 + M_2));
		ValHash += 5*(M_1 + M_2 + M_3)*(M_1 + M_2 + M_3);
		ValHash += 3*std::abs(M_1*M_2*M_3);
		
		MATRICI_SPARSE M_t(2,3);
		M_t(0,0)=L_1;
		M_t(0,1)=L_2;
		M_t(0,2)=L_3;
		M_t(1,0)=M_1;
		M_t(1,1)=M_2;
		M_t(1,2)=M_3;
		M=M_t;
		
	}
	inline unsigned long ValoreHash() const{return ValHash;}
	inline MATRICI_INTERE Mat() const {return M;}
	inline unsigned short L1() const {return L_1;}
	inline unsigned short L2() const {return L_2;}
	inline unsigned short L3() const {return L_3;}
	inline unsigned short M1() const {return M_1;}
	inline unsigned short M2() const {return M_2;}
	inline unsigned short M3() const {return M_3;}
	
};
class Hash_Numeri_3_M {
public:
	unsigned long operator()(const Numeri_3_M & N1) const {
		return N1.ValoreHash();
	}
};
class Compara_Numeri_3_M {//ordinate <H1,H2|H3>
public:
	bool operator()(const Numeri_3_M & N1 ,const Numeri_3_M & N2) {
		if ( N1.L3() == N2.L3() && N1.M3() == N2.M3()){
			if (N1.L2() == N2.L2() && N1.M2() == N2.M2()){
				if (N1.L1() == N2.L1() && N1.M1() == N2.M1()) return true;
				else return false;
			}
			else return false;
		}
		else if ( N1.L3() == N2.L3() && N1.M3() == -N2.M3()){
			if (N1.L2() == N2.L1() && N1.M2() == -N2.M1()){
				if (N1.L1() == N2.L2() && N1.M1() == -N2.M2()) return true;
				else return false;
			}
			else return false;
		}
		else return false;
	} 
};
long double Valori_3_M (const Numeri_3_M& N1){return int_3j_m(N1.L1(),N1.L2(),N1.L3(),N1.M1(),N1.M2(),N1.M3());}//ordinate <H1,H2|H3>
typedef __gnu_cxx::hash_map<Numeri_3_M, long double , Hash_Numeri_3_M, Compara_Numeri_3_M> Mappe_3_M;

class Acc_LMK { //nessuna simmetria
	unsigned short indice_N;
	std::vector <unsigned short> ll1;
	std::vector <unsigned short> ll2;
	std::vector <unsigned short> ll3;
	std::vector <short>  mm;
	unsigned long ValHash;
	MATRICI_INTERE M;
	
public:
	Acc_LMK(const MATRICI_INTERE & M_i) { //ordinate <H1,H2|H3>
		M=M_i;
		indice_N = M.size2();
		for (int i = 0; i<indice_N; ++i) {
			ll1.push_back(M(0,i));
			ll2.push_back(M(1,i));
			ll3.push_back(M(2,i));
		}

/*		ValHash = 19*(linear_sum(ll1) + linear_sum(ll2) + linear_sum(ll3));
		ValHash += 11*(linear_sum(ll1) + linear_sum(ll2) + linear_sum(ll3) )*(linear_sum(ll1) + linear_sum(ll2) + linear_sum(ll3) );
		//ValHash += 11*(linear_sum(ll1) + linear_sum(ll2) )*(linear_sum(ll1) + linear_sum(ll2) );
		//ValHash += 13*(linear_sum(ll1) + linear_sum(ll2) )*(linear_sum(ll3));
		ValHash += 7*(quad_sum(ll1) + quad_sum(ll2) + quad_sum(ll3));
		ValHash += linear_sum(ll1)*linear_sum(ll2)*linear_sum(ll3);
*/
		ValHash = 17*quad_sum(ll1);
		ValHash += 13*quad_sum(ll2);
		ValHash += 11*quad_sum(ll3);
		ValHash += linear_sum(ll1)*linear_sum(ll2)*linear_sum(ll3);
	}
	Acc_LMK(const std::vector <unsigned short> & l1,
			const std::vector <unsigned short> & l2,
			const std::vector <unsigned short> & l3){ 
		ll1 = l1;
		ll2 = l2;
		ll3 = l3;
		indice_N = ll1.size();
		MATRICI_SPARSE M_t(3,indice_N);
		for (int i = 0; i<indice_N; ++i) {
			M_t(0,i) = ll1[i];
			M_t(1,i) = ll2[i];
			M_t(2,i) = ll3[i];
		}
		M=M_t;
		
		ValHash = 17*quad_sum(ll1);
		ValHash += 13*quad_sum(ll2);
		ValHash += 11*quad_sum(ll3);
		ValHash += linear_sum(ll1)*linear_sum(ll2)*linear_sum(ll3);
		
	}

	inline unsigned long ValoreHash() const{return ValHash;}
	inline MATRICI_INTERE Mat() const {return M;}
	inline std::vector <unsigned short> l1() const {return ll1;}
	inline std::vector <unsigned short> l2() const {return ll2;}
	inline std::vector <unsigned short> l3() const {return ll3;}
	inline unsigned short I() const{return indice_N;}

};
class Hash_Acc_LMK {
public:
	unsigned long operator()(const Acc_LMK & N1) const {
		return N1.ValoreHash();
	}
};
class Compara_Acc_LMK {//ordinate <H1,H2|H3>
public:
	bool operator()(const Acc_LMK & N1 ,const Acc_LMK & N2) {
		if (N1.I() == N2.I() && delta_vec(N1.l1(),N2.l1()) && delta_vec(N1.l2(),N2.l2()) && delta_vec(N1.l3(),N2.l3()) ) return true;
		else return false;
	} 
};
typedef __gnu_cxx::hash_map<Acc_LMK, long double , Hash_Acc_LMK, Compara_Acc_LMK> Mappe_Acc_LMK;

class Integrazioni_3_P { // simple comparison function
	typedef Mappe_3_P::iterator it3P;
	typedef Mappe_3_L::iterator it3L;
	typedef Mappe_3_M::iterator it3M;
	typedef Mappe_Acc_LMK::iterator it3A;
	std::vector <Mappe_3_P> P3;
	std::vector <Mappe_3_P> P3_0;
	Mappe_3_L L3;
	Mappe_3_L L3_0;
	Mappe_3_M M3;
	Mappe_3_M M3_0;
	Mappe_Acc_LMK A_LMK;
	Mappe_Acc_LMK A_LMK_0;
	std::vector<std::string> Nomi_P3;
	std::vector<std::string> Nomi_P3_0;
	std::string Nomi_L3;
	std::string Nomi_L3_0;
	std::string Nomi_M3;
	std::string Nomi_M3_0;
	std::vector<std::string> Nomi_P3_v;
	std::vector<std::string> Nomi_P3_0_v;
	std::string Nomi_L3_v;
	std::string Nomi_L3_0_v;
	std::string Nomi_M3_v;
	std::string Nomi_M3_0_v;

	std::string Nomi_A_LMK;
	std::string Nomi_A_LMK_v;
	std::string Nomi_A_LMK_0;
	std::string Nomi_A_LMK_0_v;
	
	std::string Nome_dir;
	
	std::fstream f;
	std::fstream f_v;
	
	
public:
	
	Integrazioni_3_P(const unsigned short & NN) {
		f.precision(std::numeric_limits<long double>::digits10); 		
		f_v.precision(std::numeric_limits<long double>::digits10); 		

		P3.resize(NN - 1);
		P3_0.resize(NN - 1);
		Nome_dir = "lista_int";
		int crea_dir = mkdir(Nome_dir.c_str(), 0777);
		Nomi_L3 = "./lista_int/INT_3_L.txt";
		Nomi_L3_0 = "./lista_int/INT_3_L_0.txt";
		Nomi_M3 = "./lista_int/INT_3_M.txt";
		Nomi_M3_0 = "./lista_int/INT_3_M_0.txt";
		Nomi_A_LMK = "./lista_int/INT_A_LMK_" + scrivi(NN) + ".txt";
		Nomi_A_LMK_0 = "./lista_int/INT_A_LMK_0_" + scrivi(NN) + ".txt";
		for (short i = 1;i < NN ; ++i) {
			std::string s = "./lista_int/INT_3_P_[i_" + scrivi(i+1)+ "].txt";
			std::string s0 = "./lista_int/INT_3_P_0_[i_" + scrivi(i+1)+ "].txt";
			Nomi_P3.push_back(s);
			Nomi_P3_0.push_back(s0);			
		}		
		
		Nomi_L3_v = "./lista_int/INT_3_L_v.txt";
		Nomi_L3_0_v = "./lista_int/INT_3_L_0_v.txt";
		Nomi_M3_v = "./lista_int/INT_3_M_v.txt";
		Nomi_M3_0_v = "./lista_int/INT_3_M_0_v.txt";
		Nomi_A_LMK_v = "./lista_int/INT_A_LMK_" + scrivi(NN) + "_v.txt";
		Nomi_A_LMK_0_v = "./lista_int/INT_A_LMK_" + scrivi(NN) + "_v.txt";
		for (short i = 1;i < NN ; ++i) {
			std::string s_v = "./lista_int/INT_3_P_[i_" + scrivi(i+1)+ "]_v.txt";
			std::string s0_v = "./lista_int/INT_3_P_0_[i_" + scrivi(i+1)+ "]_v.txt";
			Nomi_P3_v.push_back(s_v);
			Nomi_P3_0_v.push_back(s0_v);			
		}		
		{//carico da file
			MATRICI_INTERE M;
			long double val; 
			for (short i = 1;i < NN ; ++i) {
				
				{					
					f.open(Nomi_P3[i-1].c_str(),std::ios::in);
					f_v.open(Nomi_P3_v[i-1].c_str(),std::ios::in);
					if (f && f_v) {
						std::cerr << "Leggo File Nomi_P3_[" << i-1 << "]" << std::endl;
						while (f >> M && f_v >> val){
							P3[i-1][Numeri_3_P(i,M)] = val;
						}
					}
					f.close();
					f_v.close();
				}
				{
					f.open(Nomi_P3_0[i-1].c_str(),std::ios::in);
					f_v.open(Nomi_P3_0_v[i-1].c_str(),std::ios::in);
					if (f && f_v) {
						std::cerr << "Leggo File Nomi_P3_0_[" << i-1 << "]" << std::endl;
						while (f >> M && f_v >> val){
							P3_0[i-1][Numeri_3_P(i,M)] = val;
						}
					}
					f.close();
					f_v.close();
				}
			}
			{
				f.open(Nomi_L3.c_str(),std::ios::in);
				f_v.open(Nomi_L3_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_L3" << std::endl;
					while (f >> M && f_v >> val){
						L3[Numeri_3_L(M)] = val;
					}
				}
				//else crea_file = true;
				f.close();
				f_v.close();
			}
			{
				f.open(Nomi_L3_0.c_str(),std::ios::in);
				f_v.open(Nomi_L3_0_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_L3_0" << std::endl;
					while (f >> M && f_v >> val){
						L3_0[Numeri_3_L(M)] = val;
					}
				}
				f.close();
				f_v.close();
			}
			{
				f.open(Nomi_M3.c_str(),std::ios::in);
				f_v.open(Nomi_M3_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_M3" << std::endl;
					while (f >> M && f_v >> val){
						M3[Numeri_3_M(M)] = val;
					}
				}
				
				f.close();
				f_v.close();
			}
			{
				f.open(Nomi_M3_0.c_str(),std::ios::in);
				f_v.open(Nomi_M3_0_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_M3_0" << std::endl;
					while (f >> M && f_v >> val){
						M3_0[Numeri_3_M(M)] = val;
					}
				}
				f.close();
				f_v.close();
			}
			{
				f.open(Nomi_A_LMK.c_str(),std::ios::in);
				f_v.open(Nomi_A_LMK_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_A_LMK" << std::endl;
					while (f >> M && f_v >> val){
						A_LMK[Acc_LMK(M)] = val;
					}
				}
				f.close();
				f_v.close();
				std::cerr << "files chiusi " << std::endl;		
			}
			{
				f.open(Nomi_A_LMK_0.c_str(),std::ios::in);
				f_v.open(Nomi_A_LMK_0_v.c_str(),std::ios::in);
				if (f && f_v) {
					std::cerr << "Leggo File Nomi_A_LMK_0" << std::endl;
					while (f >> M && f_v >> val){
						A_LMK_0[Acc_LMK(M)] = val;
					}
				}
				f.close();
				f_v.close();
			}
		}
		std::cerr << "Lista Caricata" << std::endl;

	}
	
	long double Richiesta_M(const Numeri_3_M & Rlm3){
		it3M it_0 = M3_0.find(Rlm3);
		if(it_0 == M3_0.end()){
			it3M it = M3.find(Rlm3);
			if(it != M3.end())
				return it->second;
			else {
				long double val = Valori_3_M(Rlm3);
				if (val) {
					M3[Rlm3] = val;
					f.open(Nomi_M3.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_M3_v.c_str(),std::ios::out | std::ios::app);
					f << Rlm3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				
				else {
					M3_0[Rlm3] = val;
					f.open(Nomi_M3_0.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_M3_0_v.c_str(),std::ios::out | std::ios::app);
					f << Rlm3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				return val;
			}
		}
		else return 0.;//trovato nella lista degli zeri		
		
	}
	long double Richiesta_M(const unsigned short & Indice,
							const hyper & H1,
							const hyper & H2,
							const hyper & H3){
		Numeri_3_M Rlm3(Indice,H1,H2,H3);
		it3M it_0 = M3_0.find(Rlm3);
		if(it_0 == M3_0.end()){
			it3M it = M3.find(Rlm3);
			if(it != M3.end())
				return it->second;
			else {
				long double val = Valori_3_M(Rlm3);
				if (val) {
					M3[Rlm3] = val;
					f.open(Nomi_M3.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_M3_v.c_str(),std::ios::out | std::ios::app);
					f << Rlm3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				
				else {
					M3_0[Rlm3] = 0.;
					f.open(Nomi_M3_0.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_M3_0_v.c_str(),std::ios::out | std::ios::app);
					f << Rlm3.Mat() <<std::endl;
					f_v <<  0. <<std::endl;
					f.close();
					f_v.close();
				}
				return val;
			}
		}
		else return 0.;//trovato nella lista degli zeri		
	}
	
	long double Richiesta_L(const Numeri_3_L & Rl3){
		it3L it_0 = L3_0.find(Rl3);
		if(it_0 == L3_0.end()){
			it3L it = L3.find(Rl3);
			if(it != L3.end())
				return it->second;
			else {
				long double val = Valori_3_L(Rl3);
				if (val) {
					L3[Rl3] = val;
					f.open(Nomi_L3.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_L3_v.c_str(),std::ios::out | std::ios::app);
					f << Rl3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				else {
					L3_0[Rl3] = val;
					f.open(Nomi_L3_0.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_L3_0_v.c_str(),std::ios::out | std::ios::app);
					f << Rl3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				return val;
			}
		}
		else return 0.;		
		
	}
	long double Richiesta_L(const unsigned short & Indice,
							const hyper & H1,
							const hyper & H2,
							const hyper & H3){
		Numeri_3_L Rl3(Indice,H1,H2,H3);
		
		it3L it_0 = L3_0.find(Rl3);
		if(it_0 == L3_0.end()){
			it3L it = L3.find(Rl3);
			if(it != L3.end())
				return it->second;
			else {
				long double val = Valori_3_L(Rl3);
				if (val) {
					L3[Rl3] = val;
					f.open(Nomi_L3.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_L3_v.c_str(),std::ios::out | std::ios::app);
					f << Rl3.Mat() <<std::endl;
					f_v << val <<std::endl;
					f.close();
					f_v.close();
				}
				else {
					L3_0[Rl3] =  0.;
					f.open(Nomi_L3_0.c_str(),std::ios::out | std::ios::app);
					f_v.open(Nomi_L3_0_v.c_str(),std::ios::out | std::ios::app);
					f << Rl3.Mat() <<std::endl;
					f_v <<  0. <<std::endl;
					f.close();
					f_v.close();
				}
				return val;
			}
		}
		else return 0.;		
		
	}
	
	long double Richiesta_P(const HyperInt & G, const Numeri_3_P & RP3){
		if (RP3.I()==0) return 1.;
		else{
			unsigned short II = RP3.I() - 1;
			it3P it_0 = P3_0[II].find(RP3);
			if(it_0 == P3_0[II].end()){
				it3P it = P3[II].find(RP3);
				if(it != P3[II].end())
					return it->second;
				else {
					long double val = Valori_3_P(G,RP3);
					if (std::abs(val) > errore_moltiplicazione_matrici){
						P3[II][RP3] = val;
						f.open(Nomi_P3[II].c_str(),std::ios::out | std::ios::app);
						f_v.open(Nomi_P3_v[II].c_str(),std::ios::out | std::ios::app);
						f << RP3.Mat() <<std::endl;
						f_v << val <<std::endl;
						f.close();
						f_v.close();
					}
					else {
						P3[II][RP3] = 0.;
						f.open(Nomi_P3_0[II].c_str(),std::ios::out | std::ios::app);
						f_v.open(Nomi_P3_0_v[II].c_str(),std::ios::out | std::ios::app);
						f << RP3.Mat() <<std::endl;
						f_v << 0. <<std::endl;
						f.close();
						f_v.close();
					}
					return val;
				}
			}
			else return 0.;
			
		}
	}
	long double Richiesta_P(const HyperInt & G,
							const unsigned short & Indice,
							const hyper & H1,
							const hyper & H2,
							const hyper & H3){
		if (Indice == 0) {
			return 1.;
		}
		else{
			Numeri_3_P RP3(Indice,H1,H2,H3);
			unsigned short II = RP3.I() - 1;
			it3P it_0 = P3_0[II].find(RP3);
			if(it_0 == P3_0[II].end()){
				it3P it = P3[II].find(RP3);
				if(it != P3[II].end())
					return it->second;
				else {
					long double val = Valori_3_P(G,RP3);
					if (std::abs(val) > errore_moltiplicazione_matrici){
						P3[II][RP3] = val;
						f.open(Nomi_P3[II].c_str(),std::ios::out | std::ios::app);
						f_v.open(Nomi_P3_v[II].c_str(),std::ios::out | std::ios::app);
						f << RP3.Mat() <<std::endl;
						f_v << val <<std::endl;
						f.close();
						f_v.close();
					}
					else {
						P3_0[II][RP3] = 0.;
						f.open(Nomi_P3_0[II].c_str(),std::ios::out | std::ios::app);
						f_v.open(Nomi_P3_0_v[II].c_str(),std::ios::out | std::ios::app);
						f << RP3.Mat() <<std::endl;
						f_v << 0. <<std::endl;
						f.close();
						f_v.close();
					}
					return val;
				}
			}
			else return 0.;
			
		}
	}
	
	long double Richiesta_A_LMK(bool & output,
								const std::vector <unsigned short> & l1,
								const std::vector <unsigned short> & l2,
								const std::vector <unsigned short> & l3){
		output = false;
		Acc_LMK Almk(l1,l2,l3);
		it3A it_0 = A_LMK_0.find(Almk);
		if(it_0 == A_LMK_0.end()){
			it3A it = A_LMK.find(Almk);
			if(it != A_LMK.end()){
				output = true;
				return it->second;
			}
			else {
				return 0.;
			}
		}
		else {
			output = true;
			return 0.;//trovato nella lista degli zeri		
		}
	}
	
	void Richiesta_A_LMK(const std::vector <unsigned short> & l1,
						 const std::vector <unsigned short> & l2,
						 const std::vector <unsigned short> & l3,
						 const long double & valore){
		Acc_LMK Almk(l1,l2,l3);
		if (std::abs(valore) > errore_moltiplicazione_matrici ) {
			A_LMK[Almk] = valore;
			f.open(Nomi_A_LMK.c_str(),std::ios::out | std::ios::app);
			f_v.open(Nomi_A_LMK_v.c_str(),std::ios::out | std::ios::app);
			f << Almk.Mat() <<std::endl;
			f_v << valore <<std::endl;
			f.close();
			f_v.close();
		}
		
		else {
			A_LMK_0[Almk] = 0.;
			f.open(Nomi_A_LMK_0.c_str(),std::ios::out | std::ios::app);
			f_v.open(Nomi_A_LMK_0_v.c_str(),std::ios::out | std::ios::app);
			f << Almk.Mat() <<std::endl;
			f_v << 0. <<std::endl;
			f.close();
			f_v.close();
		}
	}


};

#endif
