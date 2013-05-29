/*
 *  init_corpi.h
 *  Tesi
 *
 *  Created by Marco Resa on 10/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef INIT_CORPI_H
#define INIT_CORPI_H

#include "algebra_vettori.h"
#include "costanti_fisiche.h"

class Init_Corpi_L {
	unsigned short AA;
	unsigned short NN;
	std::vector <char> nucleoni;
	std::vector <long double> mass;
	long double mass_ref;
	std::vector <bool> charge;
	bool masse_uguali;
	
public:
	
	
	Init_Corpi_L(const std::vector <char> & nucleoni_input, const bool & masse_uguali_input){
		masse_uguali = masse_uguali_input;
		mass_ref = m_ref;
		nucleoni = nucleoni_input;
		AA = nucleoni.size();
		NN = AA-1;
		if (masse_uguali) {
			for (short i = 0; i < AA; ++i) {
				if (nucleoni[i] == 'p'){
					mass.push_back(m_ref);
					charge.push_back(1);
				}
				else if (nucleoni[i] == 'n'){
					mass.push_back(m_ref);
					charge.push_back(0);
				}
				else{
					mass.push_back(m_ref);
					charge.push_back(0);
				}
			}
		}
		else {
			for (short i = 0; i < AA; ++i) {
				if (nucleoni[i] == 'p'){
					mass.push_back(m_p);
					charge.push_back(1);
				}
				else if (nucleoni[i] == 'n'){
					mass.push_back(m_n);
					charge.push_back(0);
				}
				else{
					mass.push_back(m_ref);
					charge.push_back(0);
				}
			}
		}
	}
	
	Init_Corpi_L(const std::vector <char> & nucleoni_input){
		masse_uguali = false;
		mass_ref = m_ref;
		nucleoni = nucleoni_input;
		AA = nucleoni.size();
		NN = AA-1;
		
		for (short i = 0; i < AA; ++i) {
			if (nucleoni[i] == 'p'){
				mass.push_back(m_p);
				charge.push_back(1);
			}
			else if (nucleoni[i] == 'n'){
				mass.push_back(m_n);
				charge.push_back(0);
			}
			else{
				mass.push_back(m_ref);
				charge.push_back(0);
			}
		}
	}
	
	~Init_Corpi_L () {
	}
	
	inline unsigned short A() const {return AA;}	
	inline unsigned short N() const {return NN;}	
	inline std::vector <char> CORPI() const {return nucleoni;}	
	inline char CORPI(unsigned short i) const {return nucleoni[i];}	
	inline std::vector <long double> M() const {return mass;}	
	inline std::vector <bool> C() const {return charge;}	
	inline long double M_R() const {return mass_ref;}
	inline bool M_u() const {return masse_uguali;}	
	inline void print() const {
		std::cerr << "A = " << AA <<std::endl;
		std::cerr << "N = " << NN <<std::endl;
		std::cerr << "Corpi: " ;
		for (short i = 0; i<AA; ++i) {
			std::cerr << nucleoni[i] << " " ;			
		}
		std::cerr <<std::endl;
		std::cerr << "M_R = " << mass_ref <<std::endl;
		std::cerr << "Masse = ";
		cout_vec(mass);
		std::cerr <<std::endl;
		std::cerr << "M_u = " << masse_uguali <<std::endl;
		std::cerr << "cariche = ";
		cout_vec(charge);
		std::cerr <<std::endl;
		
	}
	inline std::string scrivi_corpi() const{
		std::string output;
		for (short i = 0; i < AA; ++i) {
			output += nucleoni[i];
		}
		return output;
	}
};

class Init_Corpi_LS {
	unsigned short AA;
	unsigned short NN;
	std::vector <char> nucleoni;
	std::vector <unsigned short> dub_s;
	short dub_s_z_min;
	std::vector <long double> mass;
	long double mass_ref;
	std::vector <bool> charge;
	bool masse_uguali;
	
public:
	
	
	Init_Corpi_LS(std::vector <char> & nucleoni_input, const bool & masse_uguali_input){
		masse_uguali = masse_uguali_input;
		mass_ref = m_ref;
		nucleoni = nucleoni_input;
		AA = nucleoni.size();
		NN = AA-1;
		if (masse_uguali) {
			for (short i = 0; i < AA; ++i) {
				if (nucleoni[i] == 'p'){
					mass.push_back(m_ref);
					charge.push_back(1);
				}
				else if (nucleoni[i] == 'n'){
					mass.push_back(m_ref);
					charge.push_back(0);
				}
				else{
					mass.push_back(m_ref);
					charge.push_back(0);
				}
				dub_s.push_back(1);
			}
		}
		else {
			for (short i = 0; i < AA; ++i) {
				if (nucleoni[i] == 'p'){
					mass.push_back(m_p);
					charge.push_back(1);
				}
				else if (nucleoni[i] == 'n'){
					mass.push_back(m_n);
					charge.push_back(0);
				}
				else{
					mass.push_back(m_ref);
					charge.push_back(0);
				}
				dub_s.push_back(1);
			}
		}
		dub_s_z_min = AA%2;
	}
	
	Init_Corpi_LS(std::vector <char> & nucleoni_input){
		masse_uguali = false;
		mass_ref = m_ref;
		nucleoni = nucleoni_input;
		AA = nucleoni.size();
		NN = AA-1;
		
		for (short i = 0; i < AA; ++i) {
			if (nucleoni[i] == 'p'){
				mass.push_back(m_p);
				charge.push_back(1);
			}
			else if (nucleoni[i] == 'n'){
				mass.push_back(m_n);
				charge.push_back(0);
			}
			else{
				mass.push_back(m_ref);
				charge.push_back(0);
			}
			dub_s.push_back(1);
		}
		dub_s_z_min = AA%2;
	}
	
	~Init_Corpi_LS () {
	}
	
	inline unsigned short A() const {return AA;}	
	inline unsigned short N() const {return NN;}	
	inline char CORPI(unsigned short i) const {return nucleoni[i];}	
	inline std::vector <long double> M() const {return mass;}	
	inline std::vector <bool> C() const {return charge;}	
	inline long double M_R() const {return mass_ref;}
	inline bool M_u() const {return masse_uguali;}
	inline std::vector <unsigned short> D_s() const{return dub_s;};
	inline short D_Z_min() const{return dub_s_z_min;};
	inline std::string scrivi_corpi() const{
		std::string output;
		for (short i = 0; i < AA; ++i) {
			output += nucleoni[i];
		}
		return output;
	}
	
};

class Init_Corpi_LST {
	unsigned short AA;
	unsigned short NN;
	std::vector <char> nucleoni;
	std::vector <vettore_ldouble> mass;
	long double mass_ref;
	std::vector <unsigned short> dub_s;
	short dub_s_z_min;
	std::vector <unsigned short> dub_t;
	short Z_tot;
	std::vector <vettore_short> dub_t_z_possibili;
	unsigned short possibilita_t_z;
	bool masse_uguali;
	
	
	
public:
	
	
	Init_Corpi_LST(const std::vector <char> & nucleoni_input, const bool & masse_uguali_input){
		masse_uguali = masse_uguali_input;
		mass_ref = m_ref;
		nucleoni = nucleoni_input;
		AA = nucleoni.size();
		NN = AA-1;
		Z_tot = 0;
		for (short i = 0; i < AA; ++i) {
			if (nucleoni[i] == 'p'){
				Z_tot+=1;
			}
			else if (nucleoni[i] == 'n'){
				Z_tot-=1;
			}
			dub_s.push_back(1);
			dub_t.push_back(1);
		}
		dub_s_z_min = AA%2;
		dub_t_z_possibili = costruisci_Z_tot(dub_s,Z_tot);
		possibilita_t_z = dub_t_z_possibili.size();
		if (masse_uguali) {
			for (short j = 0; j < possibilita_t_z; ++j) {
				vettore_ldouble mass_t;
				for (short i = 0; i < dub_t_z_possibili[j].size(); ++i) {
					mass_t.push_back(m_ref);
				}
				mass.push_back(mass_t);
			}
		}
		else {
			for (short j = 0; j < possibilita_t_z; ++j) {
				vettore_ldouble mass_t;
				for (short i = 0; i < dub_t_z_possibili[j].size(); ++i) {
					if (dub_t_z_possibili[j][i] == 1){
						mass_t.push_back(m_p);
					}
					else if (dub_t_z_possibili[j][i] == -1){
						mass_t.push_back(m_n);
					}
				}
				mass.push_back(mass_t);
			}
		}
		
	}
	
	
	~Init_Corpi_LST () {
	}
	
	inline unsigned short A() const {return AA;}	
	inline unsigned short N() const {return NN;}	
	inline std::vector <char> CORPI() const {return nucleoni;}	
	inline char CORPI(unsigned short i) const {return nucleoni[i];}	
	inline std::vector <vettore_ldouble> M() const {return mass;}	
	inline long double M_R() const {return mass_ref;}
	inline bool M_u() const {return masse_uguali;}
	inline std::vector <unsigned short> D_s() const{return dub_s;};
	inline std::vector <unsigned short> D_t() const{return dub_t;};
	inline short D_Z_min() const{return dub_s_z_min;};
	inline short D_T() const{return Z_tot;};
	inline std::vector <vettore_short> Stati_t() const{return dub_t_z_possibili;};
	inline short D() const{return possibilita_t_z;};
	inline std::string scrivi_corpi() const{
		std::string output;
		for (short i = 0; i < AA; ++i) {
			output += nucleoni[i];
		}
		return output;
	}
	
};

#endif
