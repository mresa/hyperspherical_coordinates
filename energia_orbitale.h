/*
 *  energia_orbitale.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/01/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef ENERGIA_ORBITALE_H
#define ENERGIA_ORBITALE_H

class Energia_radiale {
	
	unsigned short NN;					// numero particelle - 1
	unsigned short AA;					// numero particelle
	unsigned short n_massimo;			
	unsigned short m_massimo;
	std::vector < long double > masse;
	std::vector < bool > cariche;
	std::vector < bool > carica_coppia;
	long double massa_riferimento;
	std::vector < unsigned short > cluster;
	long double beta;
	MATRICI matrice_jacobi_inversa;
	std::vector	< vettore_ldouble > ROT_CIN;
	std::vector	< long double > C_ROT_CIN;
	std::vector < V_MATRICI_SPARSE > V_m; // vettore su n di un vettore sulle coppie
	std::vector < V_MATRICI_SPARSE > V_c_m; // vettore su n di un vettore sulle coppie
	std::vector < MATRICI_SPARSE > T_m; // vettore a 2 componenti
	std::vector < MATRICI_SPARSE > nabla_1;
	MATRICI_SPARSE R_quadro_m;
	
public:
	
	Energia_radiale (const unsigned short&  N_dato, 
					 const Gauss_Laguerre & G_m,
					 const HyperInt &  HI_m,
					 const unsigned short&  N_MAX,  
					 const unsigned short&  M_MAX,
					 const std::vector < bool > & CARICHE,
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
		cariche = CARICHE;
		for(short i = 0; i < AA; ++i)
			for(short j = i + 1; j < AA ; ++j)
				carica_coppia.push_back(cariche[i]*cariche[j]);
		
		masse = MASSE;
		massa_riferimento = MASSA_RIFERIMENTO;
		cluster = CLUSTER;
		MATRICI matrice_jacobi_diretta =  matrice_diretta(masse, massa_riferimento,cluster);
		matrice_jacobi_inversa = inversa(matrice_jacobi_diretta);
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
		}
		
		for(short i = 0; i <= n_massimo; i++){
			V_MATRICI_SPARSE V_m_t(potenziale(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_m_t);
			V_m.push_back(V_m_t); // il primo su n  il secondo sulla coppia
			V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN,carica_coppia ));
			correzione_errore(&V_c_m_t);
			V_c_m.push_back(V_c_m_t); // il primo su n  il secondo sulla coppia 
		}
		
		R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
		correzione_errore(&R_quadro_m);
	}
	
	Energia_radiale (const unsigned short&  N_dato, 
					 const Gauss_Laguerre & G_m,
					 const HyperInt &  HI_m,
					 const unsigned short&  N_MAX,  
					 const unsigned short&  M_MAX,
					 const std::vector < bool > & CARICHE,
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
		cariche = CARICHE;
		for(short i = 0; i < AA; ++i)
			for(short j = i + 1; j < AA ; ++j)
				carica_coppia.push_back(cariche[i]*cariche[j]);
		
		masse = MASSE;
		massa_riferimento = MASSA_RIFERIMENTO;
		cluster = CLUSTER;
		MATRICI matrice_jacobi_diretta =  matrice_diretta(masse, massa_riferimento,cluster);
		matrice_jacobi_inversa = inversa(matrice_jacobi_diretta);
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
		}
		
		for(short i = 0; i <= n_massimo; i++){
			V_MATRICI_SPARSE V_m_t(potenziale(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
			correzione_errore(&V_m_t);
			V_m.push_back(V_m_t); // il primo su n  il secondo sulla coppia
			V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN,carica_coppia ));
			correzione_errore(&V_c_m_t);
			V_c_m.push_back(V_c_m_t); // il primo su n  il secondo sulla coppia 
		}
		
		R_quadro_m = rho_quad_m(NN,m_massimo,G_m,beta);
		correzione_errore(&R_quadro_m);
	}
	
	~Energia_radiale () {
	}
	
	void PLUS(const unsigned short& N_plus,
			  const Gauss_Laguerre & G_m,
			  const HyperInt &  HI_m){
		if (N_plus > 0) {
			for(short i = n_massimo + 1 ; i <= n_massimo + N_plus; i++){
				
				V_MATRICI_SPARSE V_m_t(potenziale(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN ));
				correzione_errore(&V_m_t);
				V_m.push_back(V_m_t); // il primo su n  il secondo sulla coppia 
				
				V_MATRICI_SPARSE V_c_m_t(potenziale_coulombiano(NN, HI_m, G_m, m_massimo, i, beta , C_ROT_CIN, carica_coppia ));
				correzione_errore(&V_c_m_t);
				V_c_m.push_back(V_c_m_t); // il primo su n  il secondo sulla coppia 
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
	inline std::vector < bool > CARICHE() const {return cariche;}
	inline std::vector < bool > CARICHE_COPPIA() const {return carica_coppia;}
	inline long double MASSA_RIFERIMENTO() const {return massa_riferimento;}
	inline std::vector < unsigned short > CLUSTER() const {return cluster;}
	inline std::vector < vettore_ldouble > ANGOLI_ROT_CIN () const { return ROT_CIN;}
	inline std::vector < long double > COEFFICIENTE_ROT_CIN () const { return C_ROT_CIN;}
	inline std::vector < long double > ANGOLI_ROT_CIN (const unsigned short & i) const { return ROT_CIN[i];}
	inline long double COEFFICIENTE_ROT_CIN (const unsigned short & i) const { return C_ROT_CIN[i];}
	//sbagliato
/*	inline std::vector < long double > ANGOLI_ROT_CIN (const unsigned short & i, const unsigned short & j) const {
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i>j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return ROT_CIN[nn];
	}
	inline long double COEFFICIENTE_ROT_CIN ( const unsigned short &  i, const unsigned short &  j) const {
		assert (i < AA);
		assert (j < AA);
		assert (j != i);
		unsigned short ii = i;
		unsigned short jj = j;
		unsigned short nn = 0;
		if (i>j) {
			ii = j;
			jj = i;
		}
		nn = ((2*AA - ii + 1)*ii)/2 + jj - ii;
		return C_ROT_CIN[nn];
	}
 */
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
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V() const{return V_m;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_C() const{return V_c_m;} 
	inline MATRICI_SPARSE MATRICI_R_QUADRO() const{
		return R_quadro_m;
	}
	inline std::vector <MATRICI_SPARSE> NABLA() const{return nabla_1;}
};

class Energia {
	
	unsigned int degenerazione;		// numero di stati
	unsigned int degenerazione_K;	
	unsigned short n_max;	
	unsigned short n_coppie;	
	unsigned short parita;	
	std::vector <V_MATRICI_SPARSE> VV; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> VV_c; //vettore su cui sommare
	std::vector <V_MATRICI_SPARSE> TT;
	std::vector <V_MATRICI_SPARSE> EE;
	std::vector <V_MATRICI_SPARSE> EE_C;
	std::vector <V_MATRICI_SPARSE> NABLA_1;
	std::vector <Accoppiamenti_Iperangolari> TEMP_v;
	
	V_V_V_MATRICI_SPARSE R_quadro; //matrici R_i, il primo indice è di particella, il secondo su cui sommare, il terzo relativo agli spazi
	V_V_V_MATRICI_SPARSE R_delta_quadro; //matrici R_i, il primo indice è di coppia, il secondo su cui sommare, il terzo relativo agli spazi//	
	bool computo_coulomb;
	std::vector < bool > carica_coppia;
	
public:
	
	
	Energia (const HyperInt &  HI_m,
			 Integrazioni_3_P * LISTA,
			 const Energia_radiale & E_r,
			 const Accoppiamenti_Iperangolari_K_max & KK) {
		computo_coulomb=false;
		//parita = 2;
		// quello sopra è il caso più generale; vale per questo caso orbitale e diminuisce i tempi di calcolo
		// non vale nel caso di spin isospin: perché?
		parita = 2;	
		degenerazione = KK.D()*E_r.D();
		degenerazione_K = KK.D();
		n_max = E_r.N_M();
		n_coppie = E_r.MATRICI_V()[0].size();
		carica_coppia = E_r.CARICHE_COPPIA();
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
					TT1.push_back(E_r.MATRICI_T()[0]);
					TT.push_back(TT1);
					EE.push_back(TT1);
					EE_C.push_back(TT1);
				}
				{
					V_MATRICI_SPARSE TT2;
					TT2.push_back(D);
					TT2.push_back(E_r.MATRICI_T()[1]);
					TT.push_back(TT2);
					EE.push_back(TT2);
					EE_C.push_back(TT2);
				}
			}
			
			long double h_0 = 1.;
			
			{
				for(short n = 0 ; n <= E_r.N_M() ; ++n) {
					for(short i = 0 ; i < E_r.MATRICI_V()[n].size() ; ++i) {						
						std::vector <MATRICI_SPARSE> V_m_temp;
						V_m_temp.push_back(I);
						V_m_temp.push_back(h_0 * E_r.MATRICI_V()[n][i]);
						VV.push_back(V_m_temp);
						EE.push_back(V_m_temp);
						
						if(E_r.CARICHE_COPPIA()[i]){
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(I);
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C()[n][i]);
							VV_c.push_back(V_c_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_somma_temp;
							V_c_somma_temp.push_back(I);
							V_c_somma_temp.push_back(h_0 * (E_r.MATRICI_V_C()[n][i] + E_r.MATRICI_V()[n][i]));
							EE_C.push_back(V_c_somma_temp);
						}
						else {
							EE_C.push_back(V_m_temp);
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
					for (short n = 0; n<KK.N(); ++n) {
						V_MATRICI_SPARSE R_TEMP1;
						R_TEMP1.push_back(x_quad_k_t[n]);

						R_TEMP1.push_back(E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICI_R_QUADRO());
						R_TEMP.push_back(R_TEMP1);
					}
					unsigned short indice_di_coppia=0;
					for (short n = 0; n<KK.N(); ++n) {
						for (short k = n+1; k<KK.N(); ++k) {
							V_MATRICI_SPARSE R_TEMP1;
							R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
							R_TEMP1.push_back(2.*E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICE_JACOBI_INVERSA(i,k)*E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
							++indice_di_coppia;
						}
					}
					R_quadro.push_back(R_TEMP);
				}
				
				for (short i = 0; i < E_r.A(); ++i) { //su A 
					for (short j = i; j < E_r.A(); ++j) {	//sulla coppia	
						V_V_MATRICI_SPARSE R_TEMP;
						for (short n = 0; n<KK.N(); ++n) {
							V_MATRICI_SPARSE R_TEMP1;
							R_TEMP1.push_back(x_quad_k_t[n]);
							R_TEMP1.push_back((E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))*
											  (E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))
											  *E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
						}
						unsigned short indice_di_coppia=0;
						for (short n = 0; n<KK.N(); ++n) {
							for (short k = n+1; k<KK.N(); ++k) {
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
								R_TEMP1.push_back(2.*
												  (E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))*
												  (E_r.MATRICE_JACOBI_INVERSA(j,k)-E_r.MATRICE_JACOBI_INVERSA(i,k))*
												  E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
								++indice_di_coppia;
							}
						}
						R_delta_quadro.push_back(R_TEMP);
					}
				}
			}			
		}
		else if (E_r.A()>2){
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
					TT1.push_back(E_r.MATRICI_T()[0]);
					TT.push_back(TT1);
					EE.push_back(TT1);
					EE_C.push_back(TT1);
				}
				{
					V_MATRICI_SPARSE TT2;
					TT2.push_back(D);
					TT2.push_back(E_r.MATRICI_T()[1]);
					TT.push_back(TT2);
					EE.push_back(TT2);
					EE_C.push_back(TT2);
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
						NABLA_1_t.push_back(E_r.NABLA()[0]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z1);
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z2);
						NABLA_1_t.push_back(E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z3);
						NABLA_1_t.push_back(E_r.NABLA()[0] + 6 * E_r.NABLA()[1] - 12 *E_r.NABLA()[2]);
						NABLA_1.push_back(NABLA_1_t);
					}
					{
						V_MATRICI_SPARSE NABLA_1_t;
						NABLA_1_t.push_back(Z4);
						NABLA_1_t.push_back(E_r.NABLA()[1]);
						NABLA_1.push_back(NABLA_1_t);
					}		
				}
			}
			
			long double h_0 = 1./hyper0( E_r.N() - 1 );
			
			{
				for(short n = 0 ; n <= E_r.N_M() ; ++n) {
					unsigned short k_temp = 2*n;
					Accoppiamenti_Iperangolari TEMP(E_r.N(),KK.L_TOT(),KK.M(),k_temp);//,parita,E_r.CLUSTER());
					TEMP_v.push_back(TEMP);
					std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP,HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN()));
					for(short i = 0 ; i < E_r.MATRICI_V()[n].size() ; ++i) {
						MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
						correzione_errore(G_temp);
						
						std::vector <MATRICI_SPARSE> V_m_temp;
						V_m_temp.push_back(*G_temp);
						V_m_temp.push_back(h_0 * E_r.MATRICI_V()[n][i]);
						VV.push_back(V_m_temp);
						EE.push_back(V_m_temp);
						
						if(E_r.CARICHE_COPPIA()[i]){
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(*G_temp);
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C()[n][i]);
							VV_c.push_back(V_c_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_somma_temp;
							V_c_somma_temp.push_back(*G_temp);
							V_c_somma_temp.push_back(h_0 * (E_r.MATRICI_V_C()[n][i] + E_r.MATRICI_V()[n][i]));
							EE_C.push_back(V_c_somma_temp);
						}
						else {
							EE_C.push_back(V_m_temp);
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
					for (short n = 0; n<KK.N(); ++n) {
						V_MATRICI_SPARSE R_TEMP1;
						R_TEMP1.push_back(x_quad_k_t[n]);
						R_TEMP1.push_back(E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICI_R_QUADRO());
						R_TEMP.push_back(R_TEMP1);
					}
					unsigned short indice_di_coppia=0;
					for (short n = 0; n<KK.N(); ++n) {
						for (short k = n+1; k<KK.N(); ++k) {
							V_MATRICI_SPARSE R_TEMP1;
							
							R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
							R_TEMP1.push_back(2.*E_r.MATRICE_JACOBI_INVERSA(i,n)*E_r.MATRICE_JACOBI_INVERSA(i,k)*E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
							++indice_di_coppia;
						}
					}
					R_quadro.push_back(R_TEMP);
				}
				
				for (short i = 0; i < E_r.A(); ++i) { //su A 
					for (short j = i; j < E_r.A(); ++j) {	//sulla coppia	
						V_V_MATRICI_SPARSE R_TEMP;
						for (short n = 0; n<KK.N(); ++n) {
							V_MATRICI_SPARSE R_TEMP1;
							R_TEMP1.push_back(x_quad_k_t[n]);
							R_TEMP1.push_back((E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))*
											  (E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))
											  *E_r.MATRICI_R_QUADRO());
							R_TEMP.push_back(R_TEMP1);
						}
						unsigned short indice_di_coppia=0;
						for (short n = 0; n<KK.N(); ++n) {
							for (short k = n+1; k<KK.N(); ++k) {
								V_MATRICI_SPARSE R_TEMP1;
								R_TEMP1.push_back(x_ij_k_t[indice_di_coppia]);
								R_TEMP1.push_back(2.*
												  (E_r.MATRICE_JACOBI_INVERSA(j,n)-E_r.MATRICE_JACOBI_INVERSA(i,n))*
												  (E_r.MATRICE_JACOBI_INVERSA(j,k)-E_r.MATRICE_JACOBI_INVERSA(i,k))*
												  E_r.MATRICI_R_QUADRO());
								R_TEMP.push_back(R_TEMP1);
								++indice_di_coppia;
							}
						}
						R_delta_quadro.push_back(R_TEMP);
					}
				}
			}
			
		}
	}
	
	
	~Energia () {
	}
	
	
	void PLUS(const HyperInt &  HI_m,
			  Integrazioni_3_P * LISTA,
			  const Energia_radiale & E_r,
			  const Accoppiamenti_Iperangolari_K_max & KK){ //ha stesso p!!
		if (KK.D() - degenerazione_K > 0) {
			
			unsigned long numero_stati_aggiunti = KK.D() - degenerazione_K;
			std::cerr << "numero_stati_aggiunti = " << numero_stati_aggiunti << std::endl;
			{
				MATRICI_SPARSE temp_sp (accresci(TT[0][0],numero_stati_aggiunti));
				
				for (long i=degenerazione_K; i < KK.D(); ++i)
					temp_sp.push_back(i,i,1.);
				
				
				
				TT[0][0] = temp_sp;
				EE[0][0] = temp_sp;
				EE_C[0][0] = temp_sp;
				//parte scalare angolare	
				MATRICI_SPARSE temp_sp_1 (accresci(TT[1][0],numero_stati_aggiunti));
				for (long i=degenerazione_K; i < KK.D(); ++i){
					tipo_matrici val =  - KK.k(i)*(KK.k(i) + 3 * KK.N() - 2);
					if (val) temp_sp_1.push_back(i,i, val);
				} 
				
				TT[1][0]=temp_sp_1;
				EE[1][0] = temp_sp_1;
				EE_C[1][0] = temp_sp_1;
				
				
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
			{
				unsigned short indice_c = 0;
				for(short n = 0 ; n <= n_max ; ++n) {

					unsigned short k_temp = 2*n;
					std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP_v[n],HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN(),degenerazione_K));
//					Accoppiamenti_Iperangolari TEMP(E_r.N(),KK.L_TOT(),KK.M(),k_temp,parita,E_r.CLUSTER());
//					std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP,HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN(),degenerazione_K));
/*					for(short i = 0 ; i < n_coppie ; ++i) {
						
						MATRICI_SPARSE temp_sp (accresci(VV[n*n_coppie + i][0],numero_stati_aggiunti));
						
						
						MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
						correzione_errore(G_temp);
						
						temp_sp.plus_assign(*G_temp);
						
						if(E_r.CARICHE_COPPIA()[i]){
							VV_c[indice_c][0]=temp_sp;
							++indice_c;
						}
						
						VV[n*n_coppie + i][0] = temp_sp ;
						EE[n*n_coppie + i + 2][0] = temp_sp ;
						EE_C[n*n_coppie + i + 2][0] = temp_sp ;
						
					}
 */				
					for(short i = 0 ; i < n_coppie ; ++i) {
						
						MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
						correzione_errore(G_temp);
						
			
						if(E_r.CARICHE_COPPIA()[i]){
							VV_c[indice_c][0]=*G_temp;
							++indice_c;
						}
						
						VV[n*n_coppie + i][0] = *G_temp ;
						EE[n*n_coppie + i + 2][0] = *G_temp ;
						EE_C[n*n_coppie + i + 2][0] = *G_temp ;
						
					}
				}
			}

			long double h_0 = 1./hyper0( E_r.N() - 1 );
			{
				for(short n = n_max + 1 ; n <= E_r.N_M() ; ++n) {
					unsigned short k_temp = 2*n;
					Accoppiamenti_Iperangolari TEMP(E_r.N(),KK.L_TOT(),KK.M(),k_temp);//,parita,E_r.CLUSTER());
					TEMP_v.push_back(TEMP);
					std::vector<MATRICI_SPARSE> SUPERMATRICE(G_M_ij(n,TEMP,HI_m,LISTA,KK,E_r.ANGOLI_ROT_CIN()));

					for(short i = 0 ; i < E_r.MATRICI_V()[n].size() ; ++i) {

						MATRICI_SPARSE * G_temp = &(SUPERMATRICE[i]);
						correzione_errore(G_temp);
						

						std::vector <MATRICI_SPARSE> V_m_temp;
						V_m_temp.push_back(*G_temp);
						V_m_temp.push_back(h_0 * E_r.MATRICI_V()[n][i]);
						VV.push_back(V_m_temp);
						EE.push_back(V_m_temp);
						
						if(E_r.CARICHE_COPPIA()[i]){
							std::vector <MATRICI_SPARSE> V_c_m_temp;
							V_c_m_temp.push_back(*G_temp);
							V_c_m_temp.push_back(h_0 * E_r.MATRICI_V_C()[n][i]);
							VV_c.push_back(V_c_m_temp);
							
							std::vector <MATRICI_SPARSE> V_c_somma_temp;
							V_c_somma_temp.push_back(*G_temp);
							V_c_somma_temp.push_back(h_0 * (E_r.MATRICI_V_C()[n][i] + E_r.MATRICI_V()[n][i]));
							EE_C.push_back(V_c_somma_temp);
						}
						else {
							EE_C.push_back(V_m_temp);
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
			
			degenerazione = KK.D()*E_r.D();
			degenerazione_K = KK.D();
			n_max = E_r.N_M();
		}
		else;
	}
	
	inline unsigned int D() const {return degenerazione;}
	inline std::vector <V_MATRICI_SPARSE> MATRICI_T() const{return TT;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V() const{return VV;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_V_C() const{return VV_c;} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_E() const{
		return EE;
	} 
	inline std::vector <V_MATRICI_SPARSE> MATRICI_E_C() const{
		return EE_C;
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
		MATRICI V(v2m_2(x,dim_base[0]));
		MATRICI V_o(mult_2(TT[0],V));
		if (TT[1].front().nnz()) V_o.plus_assign(mult_2(TT[1],V));
		
		if (computo_coulomb) {
			for (short n = 2; n < EE_C.size(); ++n) {
				if(EE_C[n].front().nnz()){
					mult_2(EE_C[n],V,V_o);
				}
			}
		}
		
		else {
			for (short n = 2; n < EE.size(); ++n) {
				if(EE[n].front().nnz()){
					mult_2(EE[n],V,V_o);
				}
			}
		}
		mkl_free_buffers();
		return m2v(V_o);	
	}
};

namespace ietl {
	void mult(const Energia & E,
			  const VETTORI & x,
			  VETTORI & y) {
		y.assign(E.v_mult_int(x));
	}
}

#endif

