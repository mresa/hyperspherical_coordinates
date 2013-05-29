/*
 *  energia_isospinoriale.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/01/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ENERGIA_ISOSPINORIALE_H
#define ENERGIA_ISOSPINORIALE_H

#include "accoppiamenti_parita.h"
#include "pauli.h"


class Energia_radiale_isospin_masse {
	
	unsigned short NN;					// numero particelle - 1
	unsigned short AA;					// numero particelle
	unsigned short n_massimo;			
	unsigned short m_massimo;
	std::vector < long double > masse;
	long double massa_riferimento;
	std::vector < unsigned short > cluster;
	long double beta;
	MATRICI matrice_jacobi_inversa;
	std::vector	< vettore_ldouble > ROT_CIN;
	std::vector	< long double > C_ROT_CIN;
	std::vector < V_MATRICI_SPARSE > V_0; // vettore su n di un vettore sulle coppie
	std::vector < V_MATRICI_SPARSE > V_s; // vettore su n di un vettore sulle coppie
	std::vector < V_MATRICI_SPARSE > V_t; // vettore su n di un vettore sulle coppie
	std::vector < V_MATRICI_SPARSE > V_st; // vettore su n di un vettore sulle coppie
	std::vector < V_MATRICI_SPARSE > V_c; // vettore su n di un vettore sulle coppie
	std::vector < MATRICI_SPARSE > T_m; // vettore a 2 componenti
	std::vector < MATRICI_SPARSE > nabla_1;
	MATRICI_SPARSE R_quadro_m;
	
	std::vector < MATRICI_SPARSE > sigma_coppia;
	std::vector < MATRICI_SPARSE > carica_coppia;
	std::vector < MATRICI_SPARSE > tau_coppia;
	MATRICI_SPARSE ID_S;
	MATRICI_SPARSE ID_T;
	MATRICI_SPARSE tau_3_z;
	
public:
	
	Energia_radiale_isospin_masse (const Accoppiamenti_Spinoriali_S & super_base_s,
								   const Accoppiamenti_Spinoriali_S_max & super_base_t,
								   const unsigned short&  N_dato, 
								   const Gauss_Laguerre & G_m,
								   const HyperInt &  HI_m,
								   const unsigned short&  N_MAX,  
								   const unsigned short&  M_MAX, 
								   const std::vector < long double >& MASSE,  
								   const long double& MASSA_RIFERIMENTO, 
								   const std::vector < unsigned short >& CLUSTER,
								   const long double & BETA) {
		beta = BETA;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		assert (MASSE.size() == AA);
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = MASSE;
		massa_riferimento = MASSA_RIFERIMENTO;
		cluster = CLUSTER;
		matrice_jacobi_inversa = inversa(matrice_diretta(masse, massa_riferimento,cluster));
		ROT_CIN = rot_cin(matrice_jacobi_inversa);
		C_ROT_CIN = C_rot_cin(matrice_jacobi_inversa);
		
		T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
		correzione_errore(&T_m);
		
		if(AA==3){
			nabla_1.push_back(-1./2. * T_m[0]);
			nabla_1.push_back(-1./2. * T_m[1]);
			MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
			nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
			nabla_1.push_back(-1./2. * T3_temp);
			correzione_errore(&nabla_1);
			tau_3_z = tau_3_Z(super_base_t);
			correzione_errore(&tau_3_z);
		}
		
		for(short i = 0; i <= n_massimo; i++){
			
			V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_m_t[0]);
			correzione_errore(&V_m_t[1]);
			correzione_errore(&V_m_t[2]);
			correzione_errore(&V_m_t[3]);
			
			V_0.push_back(V_m_t[0]); 
			V_s.push_back(V_m_t[1]); 
			V_t.push_back(V_m_t[2]); 
			V_st.push_back(V_m_t[3]); 
			
			V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_c_m_t);
			V_c.push_back(V_c_m_t); 			
		}
		
		R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
		correzione_errore(&R_quadro_m);
		
		sigma_coppia = sigma_ij_s(super_base_s);
		correzione_errore(&sigma_coppia);
		MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
		for (long i=0; i < super_base_s.D(); ++i) {
			I.push_back(i,i,1.);
		}
		ID_S=I;
		
		carica_coppia = carica_coppia_ij(super_base_t);
		
		tau_coppia = sigma_ij_s(super_base_t);
		correzione_errore(&tau_coppia);
		MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
		for (long i=0; i < super_base_t.D(); ++i) {
			I1.push_back(i,i,1.);
		}
		ID_T=I1;
		
	}
	Energia_radiale_isospin_masse (const Accoppiamenti_Spinoriali_S & super_base_s,
								   const Accoppiamenti_Spinoriali_S_max & super_base_t,
								   const unsigned short&  N_dato, 
								   const Gauss_Laguerre & G_m,
								   const HyperInt &  HI_m,
								   const unsigned short&  N_MAX,  
								   const unsigned short&  M_MAX, 
								   const std::vector < long double >& MASSE,  
								   const long double& MASSA_RIFERIMENTO, 
								   const std::vector < unsigned short >& CLUSTER) {
		beta = beta_riferimento;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		assert (MASSE.size() == AA);
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = MASSE;
		massa_riferimento = MASSA_RIFERIMENTO;
		cluster = CLUSTER;
		matrice_jacobi_inversa = inversa(matrice_diretta(masse, massa_riferimento,cluster));
		ROT_CIN = rot_cin(matrice_jacobi_inversa);
		C_ROT_CIN = C_rot_cin(matrice_jacobi_inversa);

		T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
		correzione_errore(&T_m);
		
		if(AA==3){
			nabla_1.push_back(-1./2. * T_m[0]);
			nabla_1.push_back(-1./2. * T_m[1]);
			MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
			nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
			nabla_1.push_back(-1./2. * T3_temp);
			correzione_errore(&nabla_1);
			tau_3_z = tau_3_Z(super_base_t);
			correzione_errore(&tau_3_z);
		}
		
		for(short i = 0; i <= n_massimo; i++){
			
			V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_m_t[0]);
			correzione_errore(&V_m_t[1]);
			correzione_errore(&V_m_t[2]);
			correzione_errore(&V_m_t[3]);
			
			V_0.push_back(V_m_t[0]); 
			V_s.push_back(V_m_t[1]); 
			V_t.push_back(V_m_t[2]); 
			V_st.push_back(V_m_t[3]); 
			
			V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_c_m_t);
			V_c.push_back(V_c_m_t); 			
		}
		
		R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
		correzione_errore(&R_quadro_m);

		sigma_coppia = sigma_ij_s(super_base_s);
		correzione_errore(&sigma_coppia);
		MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
		for (long i=0; i < super_base_s.D(); ++i) {
			I.push_back(i,i,1.);
		}
		ID_S=I;

		carica_coppia = carica_coppia_ij(super_base_t);

		tau_coppia = sigma_ij_s(super_base_t);
		correzione_errore(&tau_coppia);
		MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
		for (long i=0; i < super_base_t.D(); ++i) {
			I1.push_back(i,i,1.);
		}
		ID_T=I1;

	}
	
	Energia_radiale_isospin_masse (const Accoppiamenti_Spinoriali_S_max & super_base_s,
								   const Accoppiamenti_Spinoriali_S_max & super_base_t,
								   const unsigned short&  N_dato, 
								   const Gauss_Laguerre & G_m,
								   const HyperInt &  HI_m,
								   const unsigned short&  N_MAX,  
								   const unsigned short&  M_MAX, 
								   const std::vector < long double >& MASSE,  
								   const long double& MASSA_RIFERIMENTO, 
								   const std::vector < unsigned short >& CLUSTER) {
		beta = beta_riferimento;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		assert (MASSE.size() == AA);
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = MASSE;
		massa_riferimento = MASSA_RIFERIMENTO;
		cluster = CLUSTER;
		matrice_jacobi_inversa = inversa(matrice_diretta(masse, massa_riferimento,cluster));
		ROT_CIN = rot_cin(matrice_jacobi_inversa);
		C_ROT_CIN = C_rot_cin(matrice_jacobi_inversa);
		
		T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
		correzione_errore(&T_m);
		
		if(AA==3){
			nabla_1.push_back(-1./2. * T_m[0]);
			nabla_1.push_back(-1./2. * T_m[1]);
			MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
			nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
			nabla_1.push_back(-1./2. * T3_temp);
			correzione_errore(&nabla_1);
			tau_3_z = tau_3_Z(super_base_t);
			correzione_errore(&tau_3_z);
		}
		
		for(short i = 0; i <= n_massimo; i++){
			
			V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_m_t[0]);
			correzione_errore(&V_m_t[1]);
			correzione_errore(&V_m_t[2]);
			correzione_errore(&V_m_t[3]);
			
			V_0.push_back(V_m_t[0]); 
			V_s.push_back(V_m_t[1]); 
			V_t.push_back(V_m_t[2]); 
			V_st.push_back(V_m_t[3]); 
			
			V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_c_m_t);
			V_c.push_back(V_c_m_t); 			
		}
		
		R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
		correzione_errore(&R_quadro_m);
		
		sigma_coppia = sigma_ij_s(super_base_s);
		correzione_errore(&sigma_coppia);
		MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
		for (long i=0; i < super_base_s.D(); ++i) {
			I.push_back(i,i,1.);
		}
		ID_S=I;
		
		carica_coppia = carica_coppia_ij(super_base_t);
		
		tau_coppia = sigma_ij_s(super_base_t);
		correzione_errore(&tau_coppia);
		MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
		for (long i=0; i < super_base_t.D(); ++i) {
			I1.push_back(i,i,1.);
		}
		ID_T=I1;
		
	}
	
	~Energia_radiale_isospin_masse () {
	}
	
	
	void PLUS(const unsigned short& N_plus,
			  const Gauss_Laguerre & G_m,
			  const HyperInt &  HI_m){
		if (N_plus > 0) {
			for(short i = n_massimo + 1 ; i <= n_massimo + N_plus; i++){
				V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
				correzione_errore(&V_m_t[0]);
				correzione_errore(&V_m_t[1]);
				correzione_errore(&V_m_t[2]);
				correzione_errore(&V_m_t[3]);
				V_0.push_back(V_m_t[0]); // il primo su n  il secondo sulla coppia 
				V_s.push_back(V_m_t[1]); // il primo su n  il secondo sulla coppia 
				V_t.push_back(V_m_t[2]); // il primo su n  il secondo sulla coppia 
				V_st.push_back(V_m_t[3]); // il primo su n  il secondo sulla coppia
				
				V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
				correzione_errore(&V_c_m_t);
				V_c.push_back(V_c_m_t); 			
				
			}
			n_massimo += N_plus;
		}
	}
	
	inline unsigned short A() const {return AA;}
	inline unsigned short N() const {return NN;}
	inline unsigned short N_M() const {return n_massimo;}
	inline unsigned short M_M() const {return m_massimo;}
	inline unsigned short D() const {return m_massimo + 1;}
	inline long double BETA() const {return beta;}
	inline std::vector < long double > MASSE() const {return masse;}
	inline long double MASSA_RIFERIMENTO() const {return massa_riferimento;}
	inline std::vector < unsigned short > CLUSTER() const {return cluster;}
	inline std::vector < vettore_ldouble > ANGOLI_ROT_CIN () const { return ROT_CIN;}
	inline std::vector < long double > ANGOLI_ROT_CIN (const unsigned short & i, const unsigned short & j) const {
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i<j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return ROT_CIN[nn];
	}
	inline std::vector	< long double > COEFFICIENTE_ROT_CIN() const{
		return C_ROT_CIN;
	}   
	inline long double COEFFICIENTE_ROT_CIN ( const unsigned short &  i, const unsigned short &  j) const {
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i<j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return C_ROT_CIN[nn];
	}
	inline MATRICI MATRICE_JACOBI_INVERSA() const{
		return matrice_jacobi_inversa;
	}
	inline long double MATRICE_JACOBI_INVERSA(const unsigned short &  i, const unsigned short &  j) const{
		if (j<NN) {
			return matrice_jacobi_inversa(i,AA-2-j);
		}
		else{
			return matrice_jacobi_inversa(i,NN);
		}
	}	
	inline std::vector <MATRICI_SPARSE> MATRICI_T() const{return T_m;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_0() const{return V_0;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_S() const{return V_s;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_T() const{return V_t;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_ST() const{return V_st;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_C() const{return V_c;} 
	MATRICI_SPARSE MATRICI_R_QUADRO() const{
		return R_quadro_m;
	}
	inline std::vector < MATRICI_SPARSE > SIGMA() const {return sigma_coppia;}
	inline std::vector < MATRICI_SPARSE > CARICA_COPPIA() const {return carica_coppia;}
	inline std::vector < MATRICI_SPARSE > TAU() const {return tau_coppia;}
	inline MATRICI_SPARSE I_S() const {return ID_S;}
	inline MATRICI_SPARSE I_T() const {return ID_T;}
	inline std::vector <MATRICI_SPARSE> NABLA() const{return nabla_1;}
	inline MATRICI_SPARSE TAU_3_Z() const{return tau_3_z;}
};

/*
 *  energia_isospinoriale.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/01/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

class Energia_radiale_isospin{
	
	unsigned short NN;					// numero particelle - 1
	unsigned short AA;					// numero particelle
	unsigned short n_massimo;			
	unsigned short m_massimo;
	std::vector < vettore_ldouble> masse;
	long double massa_riferimento;
	std::vector < unsigned short > cluster;
	long double beta;
	V_MATRICI matrice_jacobi_inversa;
	std::vector	< v_vettore_ldouble > ROT_CIN;	// vettore sulle rappresenzazioni
	std::vector	< vettore_ldouble > C_ROT_CIN;		// vettore sulle rappresenzazioni
	std::vector < V_V_MATRICI_SPARSE > V_0;		// vettore sulle rappresenzazioni di un vettore su n di un vettore sulle coppie
	std::vector < V_V_MATRICI_SPARSE > V_s; 
	std::vector < V_V_MATRICI_SPARSE > V_t; 
	std::vector < V_V_MATRICI_SPARSE > V_st; 
	std::vector < V_V_MATRICI_SPARSE > V_c; 
	std::vector < MATRICI_SPARSE > T_m; // vettore a 2 componenti
	std::vector < MATRICI_SPARSE > nabla_1;

	MATRICI_SPARSE R_quadro_m;
	
	std::vector < MATRICI_SPARSE > sigma_coppia;
	std::vector < MATRICI_SPARSE > carica_coppia;
	std::vector < V_MATRICI_SPARSE > tau_coppia;
	MATRICI_SPARSE ID_S;
	MATRICI_SPARSE ID_T;
	MATRICI_SPARSE tau_3_z;
	std::vector < MATRICI_SPARSE > proiettore_t;

public:
	
	
	Energia_radiale_isospin (const Accoppiamenti_Spinoriali_S & super_base_s,
							 const Accoppiamenti_Spinoriali_S_max & super_base_t,
							 const unsigned short&  N_dato, 
							 const Gauss_Laguerre & G_m,
							 const HyperInt &  HI_m,
							 const unsigned short&  N_MAX,  
							 const unsigned short&  M_MAX, 
							 const Init_Corpi_LST & Corpi,
							 const std::vector < unsigned short >& CLUSTER,
							 const long double & BETA) {
		beta = BETA;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = Corpi.M();
		massa_riferimento = Corpi.M_R();
		cluster = CLUSTER;
		
		if(Corpi.M_u()){
			Energia_radiale_isospin_masse ER_t(super_base_s,super_base_t,NN,G_m, HI_m, n_massimo, m_massimo, masse[0], massa_riferimento, cluster,beta);
			matrice_jacobi_inversa.push_back(ER_t.MATRICE_JACOBI_INVERSA());
			ROT_CIN.push_back(ER_t.ANGOLI_ROT_CIN());	// vettore sulle rappresenzazioni
			C_ROT_CIN.push_back(ER_t.COEFFICIENTE_ROT_CIN());	// vettore sulle rappresenzazioni
			V_0.push_back(ER_t.MATRICI_V_0());
			V_s.push_back(ER_t.MATRICI_V_S());
			V_t.push_back(ER_t.MATRICI_V_T());
			V_st.push_back(ER_t.MATRICI_V_ST());
			V_c.push_back(ER_t.MATRICI_V_C());
			T_m = ER_t.MATRICI_T();
			R_quadro_m = ER_t.MATRICI_R_QUADRO();
			sigma_coppia = ER_t.SIGMA();
			carica_coppia = ER_t.CARICA_COPPIA();
			tau_coppia.push_back(ER_t.TAU());
			ID_S = ER_t.I_S();
			ID_T = ER_t.I_T();
			proiettore_t.push_back(ID_T);
			if(AA==3){
				nabla_1 = ER_t.NABLA();
				tau_3_z = ER_t.TAU_3_Z();
			}
		}
		else{
			
			for (short i = 0 ; i < Corpi.D(); ++i) {
				matrice_jacobi_inversa.push_back(inversa(matrice_diretta(masse[i], massa_riferimento,cluster)));
				ROT_CIN.push_back(rot_cin(matrice_jacobi_inversa.back()));
				C_ROT_CIN.push_back(C_rot_cin(matrice_jacobi_inversa.back()));
			}
			
			T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
			correzione_errore(&T_m);
			
			if(AA==3){
				nabla_1.push_back(-1./2. * T_m[0]);
				nabla_1.push_back(-1./2. * T_m[1]);
				MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
				nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
				nabla_1.push_back(-1./2. * T3_temp);
				correzione_errore(&nabla_1);
				tau_3_z = tau_3_Z(super_base_t);
				correzione_errore(&tau_3_z);				
			}
			
			for(short casi = 0; casi < Corpi.D(); casi++){
				std::vector < V_MATRICI_SPARSE > V_0_t;
				std::vector < V_MATRICI_SPARSE > V_s_t; 
				std::vector < V_MATRICI_SPARSE > V_t_t; 
				std::vector < V_MATRICI_SPARSE > V_st_t; 
				std::vector < V_MATRICI_SPARSE > V_c_m; 
				for(short i = 0; i <= n_massimo; i++){
					
					V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_m_t[0]);
					correzione_errore(&V_m_t[1]);
					correzione_errore(&V_m_t[2]);
					correzione_errore(&V_m_t[3]);
					
					V_0_t.push_back(V_m_t[0]); 
					V_s_t.push_back(V_m_t[1]); 
					V_t_t.push_back(V_m_t[2]); 
					V_st_t.push_back(V_m_t[3]); 
					
					V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_c_m_t);
					V_c_m.push_back(V_c_m_t); 					
				}
				V_0.push_back(V_0_t); 
				V_s.push_back(V_s_t); 
				V_t.push_back(V_t_t); 
				V_st.push_back(V_st_t); 
				V_c.push_back(V_c_m);
			}
			
			R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
			correzione_errore(&R_quadro_m);
			
			sigma_coppia = sigma_ij_s(super_base_s);
			correzione_errore(&sigma_coppia);
			
			MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
			for (long i=0; i < super_base_s.D(); ++i) {
				I.push_back(i,i,1.);
			}
			ID_S = I;
			
			carica_coppia = carica_coppia_ij(super_base_t);
			
			tau_coppia = sigma_ij_s_proiettato(super_base_t,Corpi.Stati_t());
			for (long i=0; i < tau_coppia.size(); ++i) {
				correzione_errore(&(tau_coppia[i]));
			}
			
			MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
			for (long i=0; i < super_base_t.D(); ++i) {
				I1.push_back(i,i,1.);
			}
			ID_T = I1;
			
			proiettore_t = proiettore_z(super_base_t,Corpi.Stati_t());
			for (short ii = 0; ii < proiettore_t.size(); ++ii) {
				correzione_errore(&(proiettore_t[ii]));
			}
		}
	}
	Energia_radiale_isospin (const Accoppiamenti_Spinoriali_S & super_base_s,
							 const Accoppiamenti_Spinoriali_S_max & super_base_t,
							 const unsigned short&  N_dato, 
							 const Gauss_Laguerre & G_m,
							 const HyperInt &  HI_m,
							 const unsigned short&  N_MAX,  
							 const unsigned short&  M_MAX, 
							 const Init_Corpi_LST & Corpi,
							 const std::vector < unsigned short >& CLUSTER) {
		beta = beta_riferimento;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = Corpi.M();
		massa_riferimento = Corpi.M_R();
		cluster = CLUSTER;

		if(Corpi.M_u()){
			Energia_radiale_isospin_masse ER_t(super_base_s,super_base_t,NN,G_m, HI_m, n_massimo, m_massimo, masse[0], massa_riferimento, cluster);
			matrice_jacobi_inversa.push_back(ER_t.MATRICE_JACOBI_INVERSA());
			ROT_CIN.push_back(ER_t.ANGOLI_ROT_CIN());	// vettore sulle rappresenzazioni
			C_ROT_CIN.push_back(ER_t.COEFFICIENTE_ROT_CIN());	// vettore sulle rappresenzazioni
			V_0.push_back(ER_t.MATRICI_V_0());
			V_s.push_back(ER_t.MATRICI_V_S());
			V_t.push_back(ER_t.MATRICI_V_T());
			V_st.push_back(ER_t.MATRICI_V_ST());
			V_c.push_back(ER_t.MATRICI_V_C());
			T_m = ER_t.MATRICI_T();
			R_quadro_m = ER_t.MATRICI_R_QUADRO();
			sigma_coppia = ER_t.SIGMA();
			carica_coppia = ER_t.CARICA_COPPIA();
			tau_coppia.push_back(ER_t.TAU());
			ID_S = ER_t.I_S();
			ID_T = ER_t.I_T();
			proiettore_t.push_back(ID_T);
			if(AA==3){
				nabla_1 = ER_t.NABLA();
				tau_3_z = ER_t.TAU_3_Z();
			}
		}
		else{
			
			for (short i = 0 ; i < Corpi.D(); ++i) {
				matrice_jacobi_inversa.push_back(inversa(matrice_diretta(masse[i], massa_riferimento,cluster)));
				ROT_CIN.push_back(rot_cin(matrice_jacobi_inversa.back()));
				C_ROT_CIN.push_back(C_rot_cin(matrice_jacobi_inversa.back()));
			}
			
			T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
			correzione_errore(&T_m);
			
			if(AA==3){
				nabla_1.push_back(-1./2. * T_m[0]);
				nabla_1.push_back(-1./2. * T_m[1]);
				MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
				nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
				nabla_1.push_back(-1./2. * T3_temp);
				correzione_errore(&nabla_1);
				tau_3_z = tau_3_Z(super_base_t);
				correzione_errore(&tau_3_z);				
			}
			
			for(short casi = 0; casi < Corpi.D(); casi++){
				std::vector < V_MATRICI_SPARSE > V_0_t;
				std::vector < V_MATRICI_SPARSE > V_s_t; 
				std::vector < V_MATRICI_SPARSE > V_t_t; 
				std::vector < V_MATRICI_SPARSE > V_st_t; 
				std::vector < V_MATRICI_SPARSE > V_c_m; 
				for(short i = 0; i <= n_massimo; i++){
					
					V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_m_t[0]);
					correzione_errore(&V_m_t[1]);
					correzione_errore(&V_m_t[2]);
					correzione_errore(&V_m_t[3]);
					
					V_0_t.push_back(V_m_t[0]); 
					V_s_t.push_back(V_m_t[1]); 
					V_t_t.push_back(V_m_t[2]); 
					V_st_t.push_back(V_m_t[3]); 
					
					V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_c_m_t);
					V_c_m.push_back(V_c_m_t); 					
				}
				V_0.push_back(V_0_t); 
				V_s.push_back(V_s_t); 
				V_t.push_back(V_t_t); 
				V_st.push_back(V_st_t); 
				V_c.push_back(V_c_m);
			}
			
			R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
			correzione_errore(&R_quadro_m);
			
			sigma_coppia = sigma_ij_s(super_base_s);
			correzione_errore(&sigma_coppia);
			
			MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
			for (long i=0; i < super_base_s.D(); ++i) {
				I.push_back(i,i,1.);
			}
			ID_S = I;

			carica_coppia = carica_coppia_ij(super_base_t);

			tau_coppia = sigma_ij_s_proiettato(super_base_t,Corpi.Stati_t());
			for (long i=0; i < tau_coppia.size(); ++i) {
				correzione_errore(&(tau_coppia[i]));
			}

			MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
			for (long i=0; i < super_base_t.D(); ++i) {
				I1.push_back(i,i,1.);
			}
			ID_T = I1;
			
			proiettore_t = proiettore_z(super_base_t,Corpi.Stati_t());
			for (short ii = 0; ii < proiettore_t.size(); ++ii) {
				correzione_errore(&(proiettore_t[ii]));
			}
		}
	}
	Energia_radiale_isospin (const Accoppiamenti_Spinoriali_S_max & super_base_s,
							 const Accoppiamenti_Spinoriali_S_max & super_base_t,
							 const unsigned short&  N_dato, 
							 const Gauss_Laguerre & G_m,
							 const HyperInt &  HI_m,
							 const unsigned short&  N_MAX,  
							 const unsigned short&  M_MAX, 
							 const Init_Corpi_LST & Corpi,
							 const std::vector < unsigned short >& CLUSTER) {
		beta = beta_riferimento;
		assert(N_dato > 1);
		NN = N_dato;
		AA = NN + 1;
		n_massimo = N_MAX;
		m_massimo = M_MAX;
		masse = Corpi.M();
		massa_riferimento = Corpi.M_R();
		cluster = CLUSTER;
		
		if(Corpi.M_u()){
			Energia_radiale_isospin_masse ER_t(super_base_s,super_base_t,NN,G_m, HI_m, n_massimo, m_massimo, masse[0], massa_riferimento, cluster);
			matrice_jacobi_inversa.push_back(ER_t.MATRICE_JACOBI_INVERSA());
			ROT_CIN.push_back(ER_t.ANGOLI_ROT_CIN());	// vettore sulle rappresenzazioni
			C_ROT_CIN.push_back(ER_t.COEFFICIENTE_ROT_CIN());	// vettore sulle rappresenzazioni
			V_0.push_back(ER_t.MATRICI_V_0());
			V_s.push_back(ER_t.MATRICI_V_S());
			V_t.push_back(ER_t.MATRICI_V_T());
			V_st.push_back(ER_t.MATRICI_V_ST());
			V_c.push_back(ER_t.MATRICI_V_C());
			T_m = ER_t.MATRICI_T();
			R_quadro_m = ER_t.MATRICI_R_QUADRO();
			sigma_coppia = ER_t.SIGMA();
			carica_coppia = ER_t.CARICA_COPPIA();
			tau_coppia.push_back(ER_t.TAU());
			ID_S = ER_t.I_S();
			ID_T = ER_t.I_T();
			proiettore_t.push_back(ID_T);
			if(AA==3){
				nabla_1 = ER_t.NABLA();
				tau_3_z = ER_t.TAU_3_Z();
			}
		}
		else{
			
			for (short i = 0 ; i < Corpi.D(); ++i) {
				matrice_jacobi_inversa.push_back(inversa(matrice_diretta(masse[i], massa_riferimento,cluster)));
				ROT_CIN.push_back(rot_cin(matrice_jacobi_inversa.back()));
				C_ROT_CIN.push_back(C_rot_cin(matrice_jacobi_inversa.back()));
			}
			
			T_m = cinetica(NN,m_massimo,G_m,beta,massa_riferimento);
			correzione_errore(&T_m);
			
			if(AA==3){
				nabla_1.push_back(-1./2. * T_m[0]);
				nabla_1.push_back(-1./2. * T_m[1]);
				MATRICI_SPARSE T3_temp(T_m[0].size1(),T_m[0].size2());
				nabla1(NN,m_massimo,G_m,beta,massa_riferimento, &T3_temp);
				nabla_1.push_back(-1./2. * T3_temp);
				correzione_errore(&nabla_1);
				tau_3_z = tau_3_Z(super_base_t);
				correzione_errore(&tau_3_z);				
			}
			
			for(short casi = 0; casi < Corpi.D(); casi++){
				std::vector < V_MATRICI_SPARSE > V_0_t;
				std::vector < V_MATRICI_SPARSE > V_s_t; 
				std::vector < V_MATRICI_SPARSE > V_t_t; 
				std::vector < V_MATRICI_SPARSE > V_st_t; 
				std::vector < V_MATRICI_SPARSE > V_c_m; 
				for(short i = 0; i <= n_massimo; i++){
					
					V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_m_t[0]);
					correzione_errore(&V_m_t[1]);
					correzione_errore(&V_m_t[2]);
					correzione_errore(&V_m_t[3]);
					
					V_0_t.push_back(V_m_t[0]); 
					V_s_t.push_back(V_m_t[1]); 
					V_t_t.push_back(V_m_t[2]); 
					V_st_t.push_back(V_m_t[3]); 
					
					V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_c_m_t);
					V_c_m.push_back(V_c_m_t); 					
				}
				V_0.push_back(V_0_t); 
				V_s.push_back(V_s_t); 
				V_t.push_back(V_t_t); 
				V_st.push_back(V_st_t); 
				V_c.push_back(V_c_m);
			}
			
			R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
			correzione_errore(&R_quadro_m);
			
			sigma_coppia = sigma_ij_s(super_base_s);
			correzione_errore(&sigma_coppia);
			
			MATRICI_SPARSE I (super_base_s.D(),super_base_s.D());
			for (long i=0; i < super_base_s.D(); ++i) {
				I.push_back(i,i,1.);
			}
			ID_S = I;
			
			carica_coppia = carica_coppia_ij(super_base_t);
			
			tau_coppia = sigma_ij_s_proiettato(super_base_t,Corpi.Stati_t());
			for (long i=0; i < tau_coppia.size(); ++i) {
				correzione_errore(&(tau_coppia[i]));
			}
			
			MATRICI_SPARSE I1 (super_base_t.D(),super_base_t.D());
			for (long i=0; i < super_base_t.D(); ++i) {
				I1.push_back(i,i,1.);
			}
			ID_T = I1;
			
			proiettore_t = proiettore_z(super_base_t,Corpi.Stati_t());
			for (short ii = 0; ii < proiettore_t.size(); ++ii) {
				correzione_errore(&(proiettore_t[ii]));
			}
		}
	}
	

	~Energia_radiale_isospin () {
	}
	
	void PLUS(const unsigned short& N_plus,
			  const Gauss_Laguerre & G_m,
			  const HyperInt &  HI_m){
		if (N_plus > 0) {

			for(short casi = 0; casi < proiettore_t.size(); ++casi){
				for(short i = n_massimo + 1 ; i <= n_massimo + N_plus; i++){
					V_V_MATRICI_SPARSE V_m_t(potenziale_st(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] )); // primo è 0 il secondo è 1
					correzione_errore(&V_m_t[0]);
					correzione_errore(&V_m_t[1]);
					correzione_errore(&V_m_t[2]);
					correzione_errore(&V_m_t[3]);
					V_0[casi].push_back(V_m_t[0]); 
					V_s[casi].push_back(V_m_t[1]); 
					V_t[casi].push_back(V_m_t[2]); 
					V_st[casi].push_back(V_m_t[3]);
					
					V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN[casi] ));
					correzione_errore(&V_c_m_t);
					V_c[casi].push_back(V_c_m_t); 					
					
				}
			}
			n_massimo += N_plus;
		}
	}
	
	inline unsigned short A() const {return AA;}
	inline unsigned short N() const {return NN;}
	inline unsigned short N_M() const {return n_massimo;}
	inline unsigned short M_M() const {return m_massimo;}
	inline unsigned short D() const {return m_massimo + 1;}
	inline unsigned short D_T() const {return proiettore_t.size();}
	inline long double BETA() const {return beta;}
	inline std::vector < vettore_ldouble> MASSE() const {return masse;}
	inline std::vector < long double > MASSE(const unsigned short & i) const {
		assert (i < masse.size());
		return masse[i];
	}
	inline long double MASSA_RIFERIMENTO() const {return massa_riferimento;}
	inline std::vector < unsigned short > CLUSTER() const {return cluster;}
	inline std::vector	< v_vettore_ldouble > ANGOLI_ROT_CIN () const { return ROT_CIN;}
	inline std::vector	< vettore_ldouble > ANGOLI_ROT_CIN (const unsigned short & t) const { return ROT_CIN[t];}
	inline std::vector < long double > ANGOLI_ROT_CIN (const unsigned short & t, const unsigned short & i, const unsigned short & j) const {
		assert (t < ROT_CIN.size());
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i<j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return ROT_CIN[t][nn];
	}
	inline long double COEFFICIENTE_ROT_CIN (const unsigned short & t, const unsigned short &  i, const unsigned short &  j) const {
		assert (t < C_ROT_CIN.size());
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i<j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return C_ROT_CIN[t][nn];
	}
	inline V_MATRICI MATRICE_JACOBI_INVERSA() const{
		return matrice_jacobi_inversa;
	}
	inline MATRICI MATRICE_JACOBI_INVERSA(const unsigned short & t) const{
		assert(t<matrice_jacobi_inversa.size());
		return matrice_jacobi_inversa[t];
	}
	inline long double MATRICE_JACOBI_INVERSA(const unsigned short & t, const unsigned short &  i, const unsigned short &  j) const{
		assert(t<matrice_jacobi_inversa.size());
		if (j<NN) {
			return matrice_jacobi_inversa[t](i,AA-2-j);
		}
		else{
			return matrice_jacobi_inversa[t](i,NN);
		}
	}
	
	inline std::vector <MATRICI_SPARSE> MATRICI_T() const{return T_m;} 
	inline std::vector <V_V_MATRICI_SPARSE> MATRICI_V_0() const{return V_0;} 
	inline std::vector <V_V_MATRICI_SPARSE> MATRICI_V_S() const{return V_s;} 
	inline std::vector <V_V_MATRICI_SPARSE> MATRICI_V_T() const{return V_t;} 
	inline std::vector <V_V_MATRICI_SPARSE> MATRICI_V_ST() const{return V_st;} 
	inline std::vector <V_V_MATRICI_SPARSE> MATRICI_V_C() const{return V_c;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_0(const unsigned short & t) const{return V_0[t];} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_S(const unsigned short & t) const{return V_s[t];} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_T(const unsigned short & t) const{return V_t[t];} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_ST(const unsigned short & t) const{return V_st[t];} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_C(const unsigned short & t) const{return V_c[t];} 
	MATRICI_SPARSE MATRICI_R_QUADRO() const{return R_quadro_m;}
	inline std::vector < MATRICI_SPARSE > SIGMA() const {return sigma_coppia;}
	inline std::vector < MATRICI_SPARSE > CARICA_COPPIA() const {return carica_coppia;}
	inline std::vector < V_MATRICI_SPARSE > TAU() const {return tau_coppia;}
	inline std::vector < MATRICI_SPARSE > TAU(const unsigned short & t) const {return tau_coppia[t];}
	inline std::vector < MATRICI_SPARSE > P_TAU() const {return proiettore_t;}
	inline MATRICI_SPARSE P_TAU(const unsigned short & t) const {return proiettore_t[t];}
	inline MATRICI_SPARSE I_S() const {return ID_S;}
	inline MATRICI_SPARSE I_T() const {return ID_T;}
	inline std::vector <MATRICI_SPARSE> NABLA() const{return nabla_1;}
	inline MATRICI_SPARSE TAU_3_Z() const{return tau_3_z;}

};


class Energia_isospin {
	
	unsigned long degenerazione_K;	
	unsigned short parita;
	unsigned short n_max;	
	unsigned short n_coppie;	
	std::vector <V_MATRICI_SPARSE> VV_0; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> VV_S; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> VV_T; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> VV_ST; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> VV_c; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> TT;
	std::vector <V_MATRICI_SPARSE> NABLA_1;
	std::vector <Accoppiamenti_Iperangolari> TEMP_v;

	V_V_V_MATRICI_SPARSE R_quadro; //matrici R_i, il primo indice è di particella, il secondo su cui sommare, il terzo relativo agli spazi
	V_V_V_MATRICI_SPARSE R_delta_quadro; //matrici R_i, il primo indice è di coppia, il secondo su cui sommare, il terzo relativo agli spazi//	

	bool computo_coulomb;
	int fattore;

public:
		
	Energia_isospin (const HyperInt &  HI_m,
					 Integrazioni_3_P * LISTA,
					 const Energia_radiale_isospin & E_r,
					 const Accoppiamenti_Iperangolari_K_max & KK) {
		computo_coulomb=false;
		fattore = 1;
		parita = 2; //caso generale + cluster
		degenerazione_K = KK.D();
		n_max = E_r.N_M();
		n_coppie = E_r.MATRICI_V_0(0)[0].size();
		if (E_r.A() == 2) {
			MATRICI_SPARSE I (KK.D(),KK.D());
			for (long i=0; i < KK.D(); ++i) {
				I.push_back(i,i,1.);
			}
			
			{
				//parte scalare angolare	
				MATRICI_SPARSE D(KK.D(),KK.D());
				for (long i = 0; i < KK.D(); ++i) {
					int val = - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					if(val) D(i,i) = val;
				}
				{
					V_MATRICI_SPARSE TT1;
					TT1.push_back(I);
					TT1.push_back(E_r.I_S());
					TT1.push_back(E_r.I_T());				
					TT1.push_back(fattore * E_r.MATRICI_T()[0]);
					TT.push_back(TT1);
				}
				{
					V_MATRICI_SPARSE TT2;
					TT2.push_back(D);
					TT2.push_back(E_r.I_S());
					TT2.push_back(E_r.I_T());
					TT2.push_back(fattore * E_r.MATRICI_T()[1]);
					TT.push_back(TT2);
				}
				
				if(E_r.A()==3){
					MATRICI_SPARSE Z1(I);
					
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z1(i,i) *=  2 * KK.l(i)[1] * KK.l(i)[1]  -  2 * KK.l(i)[0]*KK.l(i)[0] - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					}
					correzione_errore(&Z1);
					MATRICI_SPARSE Z2(I);
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z2(i,i) *=  KK.k(i) - 2 * KK.l(i)[1]  - 2 * (KK.n(i).front()*(1 + 2*KK.n(i).front() + 2*KK.l(i)[0]) ) / (KK.k(i) + 1) ;
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 2 +  KK.k(j) - KK.n(j).front()) * (1 + KK.n(j).front()) * 2 * KK.k(j))/((2 + KK.k(j)) * (3 + KK.k(j)));
									
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()) * (4 + KK.k(j)))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z2);
					MATRICI_SPARSE Z3(Z2.size1(),Z2.size2());
					for (unsigned long i = 0; i < KK.D(); ++i) {
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*2*( 1 - KK.n(j).front()) * (2 + KK.k(j) - KK.n(j).front()) /((2 + KK.k(j)) * (3 + KK.k(j)));
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z3);
					MATRICI_SPARSE Z4(prod(D,Z3));
					correzione_errore(&Z4);
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(I);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[0]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z1);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z2);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z3);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[0] + 6 * E_r.NABLA()[1] - 12 *E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z4);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}		
				}
			}
			
			long double h_0 = fattore * 1.;
			
			{	
				for(short n = 0 ; n <= n_max ; ++n) {
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
						for(short i = 0 ; i < n_coppie ; ++i) {
							std::vector <MATRICI_SPARSE> V_0_m_temp;
							V_0_m_temp.push_back(I);
							V_0_m_temp.push_back(E_r.I_S());
							V_0_m_temp.push_back(E_r.P_TAU(tauz));
							V_0_m_temp.push_back(h_0 * E_r.MATRICI_V_0(tauz)[n][i]);
							VV_0.push_back(V_0_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(I);
							V_c_m_temp.push_back(E_r.I_S());
							V_c_m_temp.push_back(prod(E_r.P_TAU(tauz),(E_r.CARICA_COPPIA()[i])));
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C(tauz)[n][i]);
							VV_c.push_back(V_c_m_temp);
							
							std::vector <MATRICI_SPARSE> V_1_m_temp;
							V_1_m_temp.push_back(I);
							V_1_m_temp.push_back(E_r.SIGMA()[i]);
							V_1_m_temp.push_back(E_r.P_TAU(tauz));
							V_1_m_temp.push_back(h_0 * E_r.MATRICI_V_S(tauz)[n][i]);
							VV_S.push_back(V_1_m_temp);
							
							std::vector <MATRICI_SPARSE> V_2_m_temp;
							V_2_m_temp.push_back(I);
							V_2_m_temp.push_back(E_r.I_S());
							V_2_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_2_m_temp.push_back(h_0 * E_r.MATRICI_V_T(tauz)[n][i]);
							VV_T.push_back(V_2_m_temp);
							
							std::vector <MATRICI_SPARSE> V_3_m_temp;
							V_3_m_temp.push_back(I);
							V_3_m_temp.push_back(E_r.SIGMA()[i]);
							V_3_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_3_m_temp.push_back(h_0 * E_r.MATRICI_V_ST(tauz)[n][i]);
							VV_ST.push_back(V_3_m_temp);
							
						}
					}
				}
			}
			
			{
				V_MATRICI_SPARSE x_quad_k_t(x_quad_k(KK,HI_m));
				V_MATRICI_SPARSE x_ij_k_t(x_ij_k(KK,HI_m));
				correzione_errore(&x_quad_k_t);
				correzione_errore(&x_ij_k_t);
				for (short i = 0; i < E_r.A(); ++i) { //su A 
					V_V_MATRICI_SPARSE R_TEMP;
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
						for (short n = 0; n<KK.N(); ++n) {
							V_MATRICI_SPARSE R_TEMP1;
							R_TEMP1.push_back(x_quad_k_t[n]);
							R_TEMP1.push_back(E_r.I_S());
							R_TEMP1.push_back(E_r.P_TAU(tauz));
							R_TEMP1.push_back(E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
						}
						unsigned short indice_di_coppia=0;
						for (short n = 0; n<KK.N(); ++n) {
							for (short k = n+1; k<KK.N(); ++k) {
								unsigned short indice_di_coppia=0;
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
								R_TEMP1.push_back(E_r.I_S());
								R_TEMP1.push_back(E_r.P_TAU(tauz));
								R_TEMP1.push_back(2.*E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICE_JACOBI_INVERSA(tauz,i,k)*E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
								++indice_di_coppia;
							}
						}
					}
					R_quadro.push_back(R_TEMP);
					
					for (short j = i; j < E_r.A(); ++j) {	//sulla coppia	
						V_V_MATRICI_SPARSE R_TEMP;
						for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
							for (short n = 0; n<KK.N(); ++n) {
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_quad_k_t[n]);
								R_TEMP1.push_back(E_r.I_S());
								R_TEMP1.push_back(E_r.P_TAU(tauz));
								R_TEMP1.push_back((E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))*
												  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))
												  *E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
							}
							unsigned short indice_di_coppia=0;
							for (short n = 0; n<KK.N(); ++n) {
								for (short k = n+1; k<KK.N(); ++k) {
									V_MATRICI_SPARSE R_TEMP1;
									R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
									R_TEMP1.push_back(E_r.I_S());
									R_TEMP1.push_back(E_r.P_TAU(tauz));
									R_TEMP1.push_back(2.*
													  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))*
													  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,k)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,k))*
													  E_r.MATRICI_R_QUADRO());
									R_TEMP.push_back(R_TEMP1);
									++indice_di_coppia;
								}
							}
						}
						R_delta_quadro.push_back(R_TEMP);
					}
				}
			}
		}
		else if (E_r.A()>2) {
			
			{
				MATRICI_SPARSE I (KK.D(),KK.D());
				for (long i=0; i < KK.D(); ++i) {
					I.push_back(i,i,1.);
				}
				//parte scalare angolare	
				MATRICI_SPARSE D(KK.D(),KK.D());
				for (long i = 0; i < KK.D(); ++i) {
					int val = - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					if(val) D(i,i) = val;
				}
				{
					V_MATRICI_SPARSE TT1;
					TT1.push_back(I);
					TT1.push_back(E_r.I_S());
					TT1.push_back(E_r.I_T());				
					TT1.push_back(fattore * E_r.MATRICI_T()[0]);
					TT.push_back(TT1);
				}
				{
					V_MATRICI_SPARSE TT2;
					TT2.push_back(D);
					TT2.push_back(E_r.I_S());
					TT2.push_back(E_r.I_T());
					TT2.push_back(fattore * E_r.MATRICI_T()[1]);
					TT.push_back(TT2);
				}
				
				if(E_r.A()==3){
					MATRICI_SPARSE Z1(I);
					
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z1(i,i) *=  2 * KK.l(i)[1] * KK.l(i)[1]  -  2 * KK.l(i)[0]*KK.l(i)[0] - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					}
					correzione_errore(&Z1);
					MATRICI_SPARSE Z2(I);
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z2(i,i) *=  KK.k(i) - 2 * KK.l(i)[1]  - 2 * (KK.n(i).front()*(1 + 2*KK.n(i).front() + 2*KK.l(i)[0]) ) / (KK.k(i) + 1) ;
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 2 +  KK.k(j) - KK.n(j).front()) * (1 + KK.n(j).front()) * 2 * KK.k(j))/((2 + KK.k(j)) * (3 + KK.k(j)));
									
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()) * (4 + KK.k(j)))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z2);
					MATRICI_SPARSE Z3(Z2.size1(),Z2.size2());
					for (unsigned long i = 0; i < KK.D(); ++i) {
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*2*( 1 - KK.n(j).front()) * (2 + KK.k(j) - KK.n(j).front()) /((2 + KK.k(j)) * (3 + KK.k(j)));
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z3);
					MATRICI_SPARSE Z4(prod(D,Z3));
					correzione_errore(&Z4);
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(I);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[0]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z1);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z2);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z3);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[0] + 6 * E_r.NABLA()[1] - 12 *E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z4);
						NABLA_1_t.push_back(E_r.I_S());
						NABLA_1_t.push_back(E_r.TAU_3_Z());
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}		
				}
			}
			
			long double h_0 = fattore * 1./hyper0( E_r.N() - 1 );
			
			{	
				for(short n = 0 ; n <= n_max ; ++n) {
					unsigned short k_temp = 2*n;
					Accoppiamenti_Iperangolari TEMP(E_r.N(),KK.L_TOT(),KK.M(),k_temp);//,parita,E_r.CLUSTER());
					TEMP_v.push_back(TEMP);
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
						
						std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP,HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN(tauz)));
						
						for(short i = 0 ; i < n_coppie ; ++i) {
							
							MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
							correzione_errore(G_temp);
							
							
							std::vector <MATRICI_SPARSE> V_0_m_temp;
							V_0_m_temp.push_back(*G_temp);
							V_0_m_temp.push_back(E_r.I_S());
							V_0_m_temp.push_back(E_r.P_TAU(tauz));
							V_0_m_temp.push_back(h_0 * E_r.MATRICI_V_0(tauz)[n][i]);
							VV_0.push_back(V_0_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(*G_temp);
							V_c_m_temp.push_back(E_r.I_S());
							V_c_m_temp.push_back(prod(E_r.P_TAU(tauz),(E_r.CARICA_COPPIA()[i])));
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C(tauz)[n][i]);
							VV_c.push_back(V_c_m_temp);
							
							std::vector <MATRICI_SPARSE> V_1_m_temp;
							V_1_m_temp.push_back(*G_temp);
							V_1_m_temp.push_back(E_r.SIGMA()[i]);
							V_1_m_temp.push_back(E_r.P_TAU(tauz));
							V_1_m_temp.push_back(h_0 * E_r.MATRICI_V_S(tauz)[n][i]);
							VV_S.push_back(V_1_m_temp);
							
							std::vector <MATRICI_SPARSE> V_2_m_temp;
							V_2_m_temp.push_back(*G_temp);
							V_2_m_temp.push_back(E_r.I_S());
							V_2_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_2_m_temp.push_back(h_0 * E_r.MATRICI_V_T(tauz)[n][i]);
							VV_T.push_back(V_2_m_temp);
							
							std::vector <MATRICI_SPARSE> V_3_m_temp;
							V_3_m_temp.push_back(*G_temp);
							V_3_m_temp.push_back(E_r.SIGMA()[i]);
							V_3_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_3_m_temp.push_back(h_0 * E_r.MATRICI_V_ST(tauz)[n][i]);
							VV_ST.push_back(V_3_m_temp);
							
						}
					}
				}
			}
			
			{
				V_MATRICI_SPARSE x_quad_k_t(x_quad_k(KK,HI_m));
				V_MATRICI_SPARSE x_ij_k_t(x_ij_k(KK,HI_m));
				correzione_errore(&x_quad_k_t);
				correzione_errore(&x_ij_k_t);
				for (short i = 0; i < E_r.A(); ++i) { //su A 
					V_V_MATRICI_SPARSE R_TEMP;
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
						for (short n = 0; n<KK.N(); ++n) {
							V_MATRICI_SPARSE R_TEMP1;
							R_TEMP1.push_back(x_quad_k_t[n]);
							R_TEMP1.push_back(E_r.I_S());
							R_TEMP1.push_back(E_r.P_TAU(tauz));
							R_TEMP1.push_back(E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
						}
						unsigned short indice_di_coppia=0;
						for (short n = 0; n<KK.N(); ++n) {
							for (short k = n+1; k<KK.N(); ++k) {
								unsigned short indice_di_coppia=0;
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
								R_TEMP1.push_back(E_r.I_S());
								R_TEMP1.push_back(E_r.P_TAU(tauz));
								R_TEMP1.push_back(2.*E_r.MATRICE_JACOBI_INVERSA(tauz,i,n)*E_r.MATRICE_JACOBI_INVERSA(tauz,i,k)*E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
								++indice_di_coppia;
							}
						}
					}
					R_quadro.push_back(R_TEMP);
					
					for (short j = i; j < E_r.A(); ++j) {	//sulla coppia	
						V_V_MATRICI_SPARSE R_TEMP;
						for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
							for (short n = 0; n<KK.N(); ++n) {
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_quad_k_t[n]);
								R_TEMP1.push_back(E_r.I_S());
								R_TEMP1.push_back(E_r.P_TAU(tauz));
								R_TEMP1.push_back((E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))*
												  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))
												  *E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
							}
							unsigned short indice_di_coppia=0;
							for (short n = 0; n<KK.N(); ++n) {
								for (short k = n+1; k<KK.N(); ++k) {
									V_MATRICI_SPARSE R_TEMP1;
									R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
									R_TEMP1.push_back(E_r.I_S());
									R_TEMP1.push_back(E_r.P_TAU(tauz));
									R_TEMP1.push_back(2.*
													  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,n)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,n))*
													  (E_r.MATRICE_JACOBI_INVERSA(tauz,j,k)-E_r.MATRICE_JACOBI_INVERSA(tauz,i,k))*
													  E_r.MATRICI_R_QUADRO());
									R_TEMP.push_back(R_TEMP1);
									++indice_di_coppia;
								}
							}
						}
						R_delta_quadro.push_back(R_TEMP);
					}
				}
			}
		}
	}
		
	
	
	~Energia_isospin () {
	}
	
	void PLUS(const HyperInt &  HI_m,
			  Integrazioni_3_P * LISTA,
			  const Energia_radiale_isospin & E_r,
			  const Accoppiamenti_Iperangolari_K_max & KK){
		if (KK.D() - degenerazione_K > 0) {
			
			unsigned long numero_stati_aggiunti = KK.D() - degenerazione_K;
			//std::cout << "numero_stati_aggiunti = " << numero_stati_aggiunti << std::endl;
			{
				MATRICI_SPARSE temp_sp (accresci(TT[0][0],numero_stati_aggiunti));
				
				for (long i=degenerazione_K; i < KK.D(); ++i)
					temp_sp.push_back(i,i,1.);
				
				TT[0][0] = temp_sp;
				
				//parte scalare angolare	
				MATRICI_SPARSE temp_sp_1 (accresci(TT[1][0],numero_stati_aggiunti));
				for (long i=degenerazione_K; i < KK.D(); ++i){
					tipo_matrici val =  - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					if (val) temp_sp_1.push_back(i,i, val);
				} 
				
				TT[1][0]=temp_sp_1;
				
				if(E_r.A()==3){
					MATRICI_SPARSE Z1(temp_sp);
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z1(i,i) *=  2 * KK.l(i)[1] * KK.l(i)[1]  -  2 * KK.l(i)[0]*KK.l(i)[0] - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					}
					correzione_errore(&Z1);
					MATRICI_SPARSE Z2(temp_sp);
					for (unsigned long i = 0; i < KK.D(); ++i) {
						Z2(i,i) *=  KK.k(i) - 2 * KK.l(i)[1]  - 2 * (KK.n(i).front()*(1 + 2*KK.n(i).front() + 2*KK.l(i)[0]) ) / (KK.k(i) + 1) ;
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 2 +  KK.k(j) - KK.n(j).front()) * (1 + KK.n(j).front()) * 2 * KK.k(j))/((2 + KK.k(j)) * (3 + KK.k(j)));
									
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z2(i,j) = tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()) * (4 + KK.k(j)))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z2);
					MATRICI_SPARSE Z3(Z2.size1(),Z2.size2());
					for (unsigned long i = 0; i < KK.D(); ++i) {
						for(unsigned long j = 0; j < KK.D(); ++j){ // il j è il non puntato
							if(delta_vec(KK.l(i),KK.l(j)) && delta_vec(KK.L(i),KK.L(j)) ){
								if (KK.n(i).front() == (KK.n(j).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*2*( 1 - KK.n(j).front()) * (2 + KK.k(j) - KK.n(j).front()) /((2 + KK.k(j)) * (3 + KK.k(j)));
								}
								else if (KK.n(j).front() == (KK.n(i).front() + 1)) {
									Z3(i,j) = -tre_corpi_normalizzazione(KK.n(j).front(),KK.l(j)[0],KK.l(j)[1])/tre_corpi_normalizzazione(KK.n(i).front(),KK.l(i)[0],KK.l(i)[1])*(( 1 +  2*KK.l(j)[0] + 2 * KK.n(j).front()) * ( 1 +  2*KK.l(j)[1] + 2 * KK.n(j).front()))/(2 * (1 + KK.k(j)) * (2 + KK.k(j)));
								}
							}
						}
					}
					correzione_errore(&Z3);
					MATRICI_SPARSE Z4(prod(temp_sp_1,Z3));
					correzione_errore(&Z4);
					{
						NABLA_1[0][0]=temp_sp;
						NABLA_1[1][0]=Z1;
						NABLA_1[2][0]=Z2;
						NABLA_1[3][0]=Z3;
						NABLA_1[4][0]=Z4;
					}
				}
			}
			//std::cout << "numero_stati_aggiunti = " << numero_stati_aggiunti << std::endl;
			{
				unsigned short indice = 0;
				for(short n = 0 ; n <= n_max ; ++n) {
					
					unsigned short k_temp = 2*n;			
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){
						
						std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP_v[n],HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN(tauz),degenerazione_K));
						
						for(short i = 0 ; i < n_coppie ; ++i) {
							
							MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
							correzione_errore(G_temp);

							
							VV_0[indice][0] = *G_temp ;
							VV_S[indice][0] = *G_temp ;
							VV_T[indice][0] = *G_temp ;
							VV_ST[indice][0] = *G_temp ;
							VV_c[indice][0] = *G_temp ;
							++indice;
							
						}
					}
				}
			}
			long double h_0 = fattore * 1./hyper0( E_r.N() - 1 );

			{
				for(short n = n_max + 1 ; n <= E_r.N_M() ; ++n) {
					unsigned short k_temp = 2*n;
					Accoppiamenti_Iperangolari TEMP(E_r.N(),KK.L_TOT(),KK.M(),k_temp);//,parita,E_r.CLUSTER());
					TEMP_v.push_back(TEMP);
					for(short tauz = 0 ; tauz < E_r.D_T() ; ++tauz){

						std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP,HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN(tauz)));
						
						for(short i = 0 ; i < n_coppie ; ++i) {

							MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
							correzione_errore(G_temp);
						
							std::vector <MATRICI_SPARSE> V_0_m_temp;

							V_0_m_temp.push_back(*G_temp);
							V_0_m_temp.push_back(E_r.I_S());
							V_0_m_temp.push_back(E_r.P_TAU(tauz));
							V_0_m_temp.push_back(h_0 * E_r.MATRICI_V_0(tauz)[n][i]);
							VV_0.push_back(V_0_m_temp);

							std::vector <MATRICI_SPARSE> V_1_m_temp;
							V_1_m_temp.push_back(*G_temp);
							V_1_m_temp.push_back(E_r.SIGMA()[i]);
							V_1_m_temp.push_back(E_r.P_TAU(tauz));
							V_1_m_temp.push_back(h_0 * E_r.MATRICI_V_S(tauz)[n][i]);
							VV_S.push_back(V_1_m_temp);
							
							std::vector <MATRICI_SPARSE> V_2_m_temp;
							V_2_m_temp.push_back(*G_temp);
							V_2_m_temp.push_back(E_r.I_S());
							V_2_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_2_m_temp.push_back(h_0 * E_r.MATRICI_V_T(tauz)[n][i]);
							VV_T.push_back(V_2_m_temp);
							
							std::vector <MATRICI_SPARSE> V_3_m_temp;
							V_3_m_temp.push_back(*G_temp);
							V_3_m_temp.push_back(E_r.SIGMA()[i]);
							V_3_m_temp.push_back(E_r.TAU(tauz)[i]);
							V_3_m_temp.push_back(h_0 * E_r.MATRICI_V_ST(tauz)[n][i]);
							VV_ST.push_back(V_3_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(*G_temp);
							V_c_m_temp.push_back(E_r.I_S());
							V_c_m_temp.push_back(prod(E_r.P_TAU(tauz),E_r.CARICA_COPPIA()[i]));
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C(tauz)[n][i]);
							VV_c.push_back(V_c_m_temp);
							
						}
					}
				}
			}

			{
				V_MATRICI_SPARSE x_quad_k_t(x_quad_k(KK,HI_m,degenerazione_K));
				V_MATRICI_SPARSE x_ij_k_t(x_ij_k(KK,HI_m,degenerazione_K));
				correzione_errore(&x_quad_k_t);
				correzione_errore(&x_ij_k_t);

				for (short i = 0; i < R_quadro.size(); ++i) { //su A 
					for (short n = 0; n<KK.N(); ++n) {
						MATRICI_SPARSE temp_sp (accresci(R_quadro[i][n][0],numero_stati_aggiunti));
						temp_sp.plus_assign(x_quad_k_t[n]);
						R_quadro[i][n][0] = temp_sp;
					}
					
					unsigned short indice_di_coppia=0;
					for (short n = 0; n<KK.N(); ++n) {
						for (short k = n+1; k<KK.N(); ++k) {
							MATRICI_SPARSE temp_sp (accresci(R_quadro[i][KK.N()+indice_di_coppia][0],numero_stati_aggiunti));
							temp_sp.plus_assign(x_ij_k_t[indice_di_coppia]);
							R_quadro[i][KK.N()+indice_di_coppia][0] = temp_sp;								
							++indice_di_coppia;
						}
					}
				}
				
				for (short i = 0; i < R_delta_quadro.size(); ++i) { //sulla coppia
					for (short n = 0; n<KK.N(); ++n) {
						MATRICI_SPARSE temp_sp (accresci(R_delta_quadro[i][n][0],numero_stati_aggiunti));
						temp_sp.plus_assign(x_quad_k_t[n]);
						R_delta_quadro[i][n][0] = temp_sp;
					}
					unsigned short indice_di_coppia=0;
					for (short n = 0; n<KK.N(); ++n) {
						for (short k = n+1; k<KK.N(); ++k) {
							MATRICI_SPARSE temp_sp (accresci(R_delta_quadro[i][KK.N()+indice_di_coppia][0],numero_stati_aggiunti));
							temp_sp.plus_assign(x_ij_k_t[indice_di_coppia]);
							R_delta_quadro[i][KK.N()+indice_di_coppia][0] = temp_sp;
							++indice_di_coppia;
						}
					}
				}
			}
			
			degenerazione_K = KK.D();
			n_max = E_r.N_M();
		}
		else;
	}
	
	inline std::vector <V_MATRICI_SPARSE> MATRICI_T() const{return TT;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V() const{
		std::vector <V_MATRICI_SPARSE> E_temp(VV_0);
		for (short i = 0; i<VV_S.size();++i)
			E_temp.push_back(VV_S[i]);
		for (short i = 0; i<VV_T.size();++i)
			E_temp.push_back(VV_T[i]);
		for (short i = 0; i<VV_ST.size();++i)
			E_temp.push_back(VV_ST[i]);
		return E_temp;
	}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V0() const{
		return VV_0;
	}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V1() const{
		return VV_S;
	}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V2() const{
		return VV_T;
	}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V3() const{
		return VV_ST;
	}
	
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_C() const{
		return VV_c;
	}
 	inline std::vector <V_MATRICI_SPARSE> MATRICI_E() const{
		std::vector <V_MATRICI_SPARSE> E_temp(TT);
		for (short i = 0; i<VV_0.size();++i)
			E_temp.push_back(VV_0[i]);
		for (short i = 0; i<VV_S.size();++i)
			E_temp.push_back(VV_S[i]);
		for (short i = 0; i<VV_T.size();++i)
			E_temp.push_back(VV_T[i]);
		for (short i = 0; i<VV_ST.size();++i)
			E_temp.push_back(VV_ST[i]);
		return E_temp;
	} 
 	inline std::vector <V_MATRICI_SPARSE> MATRICI_E_C() const{
		std::vector <V_MATRICI_SPARSE> E_temp(TT);
		for (short i = 0; i<VV_0.size();++i)
			E_temp.push_back(VV_0[i]);
		for (short i = 0; i<VV_S.size();++i)
			E_temp.push_back(VV_S[i]);
		for (short i = 0; i<VV_T.size();++i)
			E_temp.push_back(VV_T[i]);
		for (short i = 0; i<VV_ST.size();++i)
			E_temp.push_back(VV_ST[i]);
		for (short i = 0; i<VV_c.size();++i)
			E_temp.push_back(VV_c[i]);
		return E_temp;
	}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_NABLA() const{return NABLA_1;} 
	V_V_V_MATRICI_SPARSE MATRICI_R_QUADRO() const{
		return R_quadro;
	}
	V_V_V_MATRICI_SPARSE MATRICI_R_DELTA_QUADRO() const{
		return R_delta_quadro;
	}
	void COMPUTA_COULOMB(){
		computo_coulomb = true;
	}
	void N_COMPUTA_COULOMB(){
		computo_coulomb = false;
	}
	VETTORI v_mult_int(const VETTORI & x) const{
		vettore_ushort dim_base;
		for (short i = 0; i<TT[0].size(); ++i) {
			dim_base.push_back(TT[0][i].size1());
		}
		
		V_V_MATRICI V(v2m_4(x,dim_base));
		MATRICI V_o(mult_4(TT[0],V));
		if (TT[1].front().nnz()) V_o.plus_assign(mult_4(TT[1],V));
		
		if(computo_coulomb){
			
			for (short n = 0; n < VV_0.size(); ++n) { //di somma
				if(VV_0[n].front().nnz()){
					
					MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]*dim_base[3]);
#pragma omp parallel for schedule( guided )
					for ( long i = 0; i < dim_base[0]; ++i) {
						MATRICI M1((mult_3(VV_0[n][1],VV_0[n][2],VV_0[n][3],V[i]) + mult_3(VV_S[n][1],VV_S[n][2],VV_S[n][3],V[i]) + mult_3(VV_T[n][1],VV_T[n][2],VV_T[n][3],V[i]) + mult_3(VV_ST[n][1],VV_ST[n][2],VV_ST[n][3],V[i]) + mult_3(VV_c[n][1],VV_c[n][2],VV_c[n][3],V[i])));
						row(M2, i) = m2v(M1);
					}
					
					prod_mkl_sparse(VV_0[n].front(),M2,V_o);
				}
			}
		}
		
		else{
			
			for (short n = 0; n < VV_0.size(); ++n) { //di somma
				if(VV_0[n].front().nnz()){
					
					MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]*dim_base[3]);
#pragma omp parallel for schedule( guided )
					for ( long i = 0; i < dim_base[0]; ++i) {
						MATRICI M1((mult_3(VV_0[n][1],VV_0[n][2],VV_0[n][3],V[i]) + mult_3(VV_S[n][1],VV_S[n][2],VV_S[n][3],V[i]) + mult_3(VV_T[n][1],VV_T[n][2],VV_T[n][3],V[i]) + mult_3(VV_ST[n][1],VV_ST[n][2],VV_ST[n][3],V[i])));
						row(M2, i) = m2v(M1);
					}
					
					prod_mkl_sparse(VV_0[n].front(),M2,V_o);
				}
			}
		}
		
		mkl_free_buffers();
		return m2v(V_o);	
	}
};

namespace ietl {
	void mult(const Energia_isospin & E,
			  const VETTORI & x,
			  VETTORI & y) {
		y.assign(E.v_mult_int(x));
	}
}


#endif


