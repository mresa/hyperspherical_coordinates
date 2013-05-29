/*
 *  pauli.h
 *  Tesi
 *	VALIDO TRA STATI 1/2
 *  Created by Marco Resa on 27/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef PAULI_H
#define PAULI_H
#include <vector>
#include "algebra_matrici.h"
#include "accoppiamenti_spinoriali.h"
#include <omp.h> //calcolo parallelo OMP



//Cosa fa? matrici sigma (x,y,z) tra stati dub_z1 e dub_z2, a 3 componenti
inline std::vector <short> sigma_s (const short & dub_z1, const short & dub_z2){
	assert(dub_z1 == 1 || dub_z1 == -1 );
	assert(dub_z2 == 1 || dub_z2 == -1 );
	std::vector <short> sigma_s_t;
	
	if (dub_z1 != dub_z2){
		sigma_s_t.push_back(1);
		sigma_s_t.push_back(dub_z2); //attenzione!! c'è una i a moltiplicare
		sigma_s_t.push_back(0);
	}
	else {
		sigma_s_t.push_back(0);
		sigma_s_t.push_back(0);
		sigma_s_t.push_back(dub_z2);
	}
	return sigma_s_t;
}

//prodotto scalare sigma_i . sigma_j tra stati nella base z_1 … z_A
//è simmetrica sui vettori di base ? 
//	sì poiché è reale per costruzione

inline short sigma_ij_s (const std::vector<short> & dub_z1,
						 const std::vector<short> & dub_z2,
						 const unsigned short & indice,
						 const unsigned short & jndice){
	assert(indice != jndice);
	assert(dub_z1.size() == dub_z2.size());
	bool verita = true;
	for (short i = 0; i < dub_z1.size(); ++i) {
		assert(dub_z1[i] == 1 || dub_z1[i] == -1);
		assert(dub_z2[i] == 1 || dub_z2[i] == -1);
		if( i != indice && i != jndice) verita *= (dub_z1[i] == dub_z2[i]);
	}
	short risultato = 0.;
	if (verita){
		std::vector <short> sigma_s_i(sigma_s(dub_z1[indice],dub_z2[indice]));
		std::vector <short> sigma_s_j(sigma_s(dub_z1[jndice],dub_z2[jndice]));
		risultato += sigma_s_i[0]*sigma_s_j[0];
		risultato -= sigma_s_i[1]*sigma_s_j[1]; //dovuto alla i !!!
		risultato += sigma_s_i[2]*sigma_s_j[2];
	}
	return risultato;
}


inline V_MATRICI_SPARSE sigma_ij_s(const Accoppiamenti_Spinoriali_S_max & TEMP){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	std::vector<MATRICI_SPARSE> SUPERMATRICE;
	
	for (short indice = 0; indice < A ; ++indice) {
		for (short jndice = indice + 1; jndice < A ; ++jndice) {
			MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
			
			for (short r = 0; r < MATRI.size1(); ++r) {
				for (short c = r; c < MATRI.size2(); ++c) {
					long double C = 0.;
					for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
						long double C1 = 0.;
						for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
							C1 += TEMP.cm(c,c1) * sigma_ij_s(TEMP.D_z(r,r1),TEMP.D_z(c,c1),indice,jndice);
						}
						C += C1*TEMP.cm(r,r1);						
					}
					if(std::abs(C) > errore_matrici){
						MATRI(r,c) = C;  
						MATRI(c,r) = C;
					}
				}
			}
			SUPERMATRICE.push_back(MATRI);
		}
	}
	unsigned short S = ((A-1)*A)/2;
	//	verifico che la dimensione sia pari al numero delle coppie
	assert(SUPERMATRICE.size() == S);
	return SUPERMATRICE;
}

inline V_MATRICI_SPARSE sigma_ij_s(const Accoppiamenti_Spinoriali_S & TEMP){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	std::vector<MATRICI_SPARSE> SUPERMATRICE;
	
	for (short indice = 0; indice < A ; ++indice) {
		for (short jndice = indice + 1; jndice < A ; ++jndice) {
			MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
			
			for (short r = 0; r < MATRI.size1(); ++r) {
				for (short c = r; c < MATRI.size2(); ++c) {
					long double C = 0.;
					for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
						long double C1 = 0.;
						for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
							C1 += TEMP.cm(c,c1) * sigma_ij_s(TEMP.D_z(r,r1),TEMP.D_z(c,c1),indice,jndice);
						}
						C += C1*TEMP.cm(r,r1);						
					}
					if(std::abs(C) > errore_matrici){
						MATRI(r,c) = C;  
						MATRI(c,r) = C;
					}
				}
			}
			SUPERMATRICE.push_back(MATRI);
		}
	}
	unsigned short S = ((A-1)*A)/2;
	assert(SUPERMATRICE.size() == S);
	return SUPERMATRICE;
}

inline MATRICI_SPARSE proiettore_z(const Accoppiamenti_Spinoriali_S & TEMP, 
								   const std::vector < short > & dub_zz){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
	
	for (short r = 0; r < MATRI.size1(); ++r) {
		for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
			if (delta_vec(TEMP.D_z(r,r1), dub_zz)) {
				for (short c = r; c < MATRI.size2(); ++c) {
					for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
						if (delta_vec(TEMP.D_z(c,c1), dub_zz)) {
							long double C = TEMP.cm(r,r1)*TEMP.cm(c,c1);
							if(std::abs(C) > errore_matrici){
								MATRI(r,c) = C; 
								MATRI(c,r) = C;
							}
						}
					}
				}
			}
		}
	}
	return MATRI;
}
//indicizzato sulle possibilità
inline V_MATRICI_SPARSE proiettore_z(const Accoppiamenti_Spinoriali_S & TEMP,
									 const std::vector <vettore_short> & possibilita){
	V_MATRICI_SPARSE proiettore_z_t;
	for (short i = 0; i < possibilita.size(); ++i) {
		proiettore_z_t.push_back(proiettore_z(TEMP,possibilita[i]));
	}
	return proiettore_z_t;
}


//sulle coppie
inline V_MATRICI_SPARSE sigma_ij_s_proiettato(const Accoppiamenti_Spinoriali_S & TEMP, 
											  const std::vector < short > & dub_zz){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	std::vector<MATRICI_SPARSE> SUPERMATRICE;
	
	for (short indice = 0; indice < A ; ++indice) {
		for (short jndice = indice + 1; jndice < A ; ++jndice) {
			MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
			
			for (short r = 0; r < MATRI.size1(); ++r) {
				for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
					if (delta_vec(TEMP.D_z(r,r1), dub_zz)) {
						long double C = TEMP.cm(r,r1);
						for (short c = r; c < MATRI.size2(); ++c) {
							long double C1 = 0;
							for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
								C1 += TEMP.cm(c,c1) * sigma_ij_s(dub_zz,TEMP.D_z(c,c1),indice,jndice);
							}
							C1*=C;
							if(std::abs(C1) > errore_matrici){
								MATRI(r,c) = C1; 
								MATRI(c,r) = C1;
							}
						}
					}
				}
			}
			SUPERMATRICE.push_back(MATRI);
		}
	}
	unsigned short S = ((A-1)*A)/2;
	assert(SUPERMATRICE.size() == S);
	return SUPERMATRICE;
}

//indicizzato sulle possibilità
inline V_V_MATRICI_SPARSE sigma_ij_s_proiettato(const Accoppiamenti_Spinoriali_S & TEMP,
												const std::vector <vettore_short> & possibilita){
	V_V_MATRICI_SPARSE sigma_ij_s_proiettato_t;
	for (short i = 0; i < possibilita.size(); ++i) {
		sigma_ij_s_proiettato_t.push_back(sigma_ij_s_proiettato(TEMP,possibilita[i]));
	}
	return sigma_ij_s_proiettato_t;
}

inline MATRICI_SPARSE proiettore_z(const Accoppiamenti_Spinoriali_S_max & TEMP, 
								   const std::vector < short > & dub_zz){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
	
	for (short r = 0; r < MATRI.size1(); ++r) {
		for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
			if (delta_vec(TEMP.D_z(r,r1), dub_zz)) {
				for (short c = r; c < MATRI.size2(); ++c) {
					for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
						if (delta_vec(TEMP.D_z(c,c1), dub_zz)) {
							long double C = TEMP.cm(r,r1)*TEMP.cm(c,c1);
							if(std::abs(C) > errore_matrici){
								MATRI(r,c) = C; 
								MATRI(c,r) = C;
							}
						}
					}
				}
			}
		}
	}
	return MATRI;
}
//indicizzato sulle possibilità
inline V_MATRICI_SPARSE proiettore_z(const Accoppiamenti_Spinoriali_S_max & TEMP,
									 const std::vector <vettore_short> & possibilita){
	V_MATRICI_SPARSE proiettore_z_t;
	for (short i = 0; i < possibilita.size(); ++i) {
		proiettore_z_t.push_back(proiettore_z(TEMP,possibilita[i]));
		std::cerr << proiettore_z(TEMP,possibilita[i]) << std::endl;
	}
	return proiettore_z_t;
}


//sulle coppie
inline V_MATRICI_SPARSE sigma_ij_s_proiettato(const Accoppiamenti_Spinoriali_S_max & TEMP, 
											  const std::vector < short > & dub_zz){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	std::vector<MATRICI_SPARSE> SUPERMATRICE;
	
	for (short indice = 0; indice < A ; ++indice) {
		for (short jndice = indice + 1; jndice < A ; ++jndice) {
			MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
			
			for (short r = 0; r < MATRI.size1(); ++r) {
				for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
					if (delta_vec(TEMP.D_z(r,r1), dub_zz)) {
						long double C = TEMP.cm(r,r1);
						for (short c = r; c < MATRI.size2(); ++c) {
							long double C1 = 0;
							for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
								C1 += TEMP.cm(c,c1) * sigma_ij_s(dub_zz,TEMP.D_z(c,c1),indice,jndice);
							}
							C1*=C;
							if(std::abs(C1) > errore_matrici){
								MATRI(r,c) = C1; 
								MATRI(c,r) = C1;
							}
						}
					}
				}
			}
			SUPERMATRICE.push_back(MATRI);
		}
	}
	unsigned short S = ((A-1)*A)/2;
	assert(SUPERMATRICE.size() == S);
	return SUPERMATRICE;
}

//indicizzato sulle possibilità
inline V_V_MATRICI_SPARSE sigma_ij_s_proiettato(const Accoppiamenti_Spinoriali_S_max & TEMP,
												const std::vector <vettore_short> & possibilita){
	V_V_MATRICI_SPARSE sigma_ij_s_proiettato_t;
	for (short i = 0; i < possibilita.size(); ++i) {
		sigma_ij_s_proiettato_t.push_back(sigma_ij_s_proiettato(TEMP,possibilita[i]));
	}
	return sigma_ij_s_proiettato_t;
}

inline short carica_elettrica (const short & dub_t1, const short & dub_t2){
	assert(dub_t1 == 1 || dub_t1 == -1 );
	assert(dub_t2 == 1 || dub_t2 == -1 );
	short c_e = 0;	
	if (dub_t2 == 1 && dub_t2 == 1){
		c_e = 1;
	}
	return c_e;
}

inline short carica_coppia_ij (const std::vector<short> & dub_z1,
							   const std::vector<short> & dub_z2,
							   const unsigned short & indice,
							   const unsigned short & jndice){
	assert(indice != jndice);
	assert(dub_z1.size() == dub_z2.size());
	bool verita = true;
	for (short i = 0; i < dub_z1.size(); ++i) {
		assert(dub_z1[i] == 1 || dub_z1[i] == -1);
		assert(dub_z2[i] == 1 || dub_z2[i] == -1);
		if( i != indice && i != jndice) verita *= (dub_z1[i] == dub_z2[i]);
	}
	short risultato = 0.;
	if (verita){
		risultato = carica_elettrica(dub_z1[indice],dub_z2[indice])*carica_elettrica(dub_z1[jndice],dub_z2[jndice]);
	}
	return risultato;
}

inline V_MATRICI_SPARSE carica_coppia_ij(const Accoppiamenti_Spinoriali_S_max & TEMP){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	//costruisco un vettore di matrici, ogni elemento si riferisce alla coppia ij considerata nel mio ordinamento
	
	std::vector<MATRICI_SPARSE> SUPERMATRICE;
	
	for (short indice = 0; indice < A ; ++indice) {
		for (short jndice = indice + 1; jndice < A ; ++jndice) {
			MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
			
			for (short r = 0; r < MATRI.size1(); ++r) {
				for (short c = r; c < MATRI.size2(); ++c) {
					long double C = 0.;
					for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
						long double C1 = 0.;
						for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
							C1 += TEMP.cm(c,c1) * carica_coppia_ij(TEMP.D_z(r,r1),TEMP.D_z(c,c1),indice,jndice);
						}
						C += C1*TEMP.cm(r,r1);						
					}
					if(std::abs(C) > errore_matrici){
						MATRI(r,c) = C;  
						MATRI(c,r) = C;
					}
				}
			}
			SUPERMATRICE.push_back(MATRI);
		}
	}
	unsigned short S = ((A-1)*A)/2;
	//	verifico che la dimensione sia pari al numero delle coppie
	assert(SUPERMATRICE.size() == S);
	return SUPERMATRICE;
}

inline MATRICI_SPARSE tau_3_Z(const Accoppiamenti_Spinoriali_S & TEMP){
	unsigned short A = TEMP.A();	//Numero di particelle
		
	
	MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
	
	if (A==3) {
		for (short r = 0; r < MATRI.size1(); ++r) {
			for (short c = r; c < MATRI.size2(); ++c) {
				long double C = 0.;
				for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
					for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
						if (delta_vec(TEMP.D_z(r,r1),TEMP.D_z(c,c1))) {
							C += TEMP.cm(r,r1)*TEMP.cm(c,c1)*TEMP.D_z(c,c1)[2];
						}
					}
				}
				if(std::abs(C) > errore_matrici){
					MATRI(r,c) = C;  
					MATRI(c,r) = C; 
				}
				
			}
		}
	}
	return MATRI;
}

inline MATRICI_SPARSE tau_3_Z(const Accoppiamenti_Spinoriali_S_max & TEMP){
	unsigned short A = TEMP.A();	//Numero di particelle
	
	MATRICI_SPARSE MATRI (TEMP.D(),TEMP.D());
	
	if (A==3) {
		for (short r = 0; r < MATRI.size1(); ++r) {
			for (short c = r; c < MATRI.size2(); ++c) {
				long double C = 0.;
				for (short r1 = 0; r1 < TEMP.Dz(r); ++r1) {
					for (short c1 = 0; c1 < TEMP.Dz(c); ++c1) {
						if (delta_vec(TEMP.D_z(r,r1),TEMP.D_z(c,c1))) {
							C += TEMP.cm(r,r1)*TEMP.cm(c,c1)*TEMP.D_z(c,c1)[2];
						}
					}
				}
				if(std::abs(C) > errore_matrici){
					MATRI(r,c) = C;  
					MATRI(c,r) = C; 
				}
				
			}
		}
	}
	return MATRI;
}


#endif
