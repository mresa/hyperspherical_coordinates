/*
 *  accoppiamenti_parita.h
 *  Tesi
 *
 *  Created by Marco Resa on 01/02/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#ifndef ACC_PARITA_H
#define ACC_PARITA_H

#include "algebra_matrici.h"
#include "accoppiamenti_iperangolari.h"
#include "accoppiamenti_spinoriali.h"

class Accoppiamenti_Parita {

	typedef MATRICI_SPARSE::const_iterator1 i1_t;
	typedef MATRICI_SPARSE::const_iterator2 i2_t;
	unsigned long dimensione_base_totale;
	unsigned long dimensione_base_pari;
	unsigned long dimensione_base_dispari;
	unsigned long dimensione_base_iperangolare;
	unsigned long dimensione_base_iperangolare_pari;
	unsigned long dimensione_base_iperangolare_dispari;
	std::vector < unsigned short > parita_della_base_totale;							
	std::vector < unsigned short > parita_della_base_iperangolare;							
	std::vector < unsigned long > indice_p;							
	std::vector < unsigned long > indice_d;							
	
public:

	Accoppiamenti_Parita(const Accoppiamenti_Iperangolari_K_max & KK, 
						 const unsigned short & dimensione_m){
		dimensione_base_pari = 0;
		dimensione_base_iperangolare_pari = 0;
		dimensione_base_iperangolare_dispari = 0;
		dimensione_base_iperangolare = KK.D();
		dimensione_base_totale = dimensione_base_iperangolare*dimensione_m;
		unsigned short n_cluster = KK.P_l(0).size();
		for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
			unsigned short pp = linear_sum_bool(KK.P_l(k_i));
			if(pp == 0){//tutti pari
				parita_della_base_iperangolare.push_back(0);
				dimensione_base_iperangolare_pari++;
				indice_p.push_back(dimensione_base_iperangolare_pari);
				indice_d.push_back(0);
			}
			else if (pp == n_cluster){//tutti dispari
				parita_della_base_iperangolare.push_back(1);
				dimensione_base_iperangolare_dispari++;
				indice_p.push_back(0);
				indice_d.push_back(dimensione_base_iperangolare_dispari);
			}
			else{//misto
				parita_della_base_iperangolare.push_back(2);
				indice_p.push_back(0);
				indice_d.push_back(0);
			}
			for (short m_i = 0 ; m_i<dimensione_m; ++m_i) {
				parita_della_base_totale.push_back(parita_della_base_iperangolare.back());
			}
		}
		dimensione_base_pari = dimensione_base_iperangolare_pari * dimensione_m;
		dimensione_base_dispari = dimensione_base_iperangolare_dispari * dimensione_m;
	}
	
	Accoppiamenti_Parita(const Accoppiamenti_Iperangolari_K_max & KK, 
						 const Accoppiamenti_Spinoriali_S & SS, 
						 const unsigned short & dimensione_m){
		dimensione_base_pari = 0;
		dimensione_base_iperangolare_pari = 0;
		dimensione_base_iperangolare_dispari = 0;
		dimensione_base_iperangolare = KK.D()*SS.D();
		dimensione_base_totale = dimensione_base_iperangolare*dimensione_m;
		unsigned short n_cluster = KK.P_l()[0].size();
		for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
			for (short s_i = 0 ; s_i<SS.D(); ++s_i) {
				std::vector<unsigned short> p_temp;
				for (short i = 0 ; i<n_cluster; ++i) {
					p_temp.push_back((SS.P_S()[s_i][i] + KK.P_l()[k_i][i])%2);
				}
				unsigned short pp = linear_sum(p_temp);
				if(pp == 0){//tutti pari
					parita_della_base_iperangolare.push_back(0);
					dimensione_base_iperangolare_pari++;
					indice_p.push_back(dimensione_base_iperangolare_pari);
					indice_d.push_back(0);
				}
				else if (pp == n_cluster){//tutti dispari
					parita_della_base_iperangolare.push_back(1);
					dimensione_base_iperangolare_dispari++;
					indice_p.push_back(0);
					indice_d.push_back(dimensione_base_iperangolare_dispari);
				}
				else{//misto
					parita_della_base_iperangolare.push_back(2);
					indice_p.push_back(0);
					indice_d.push_back(dimensione_base_iperangolare_dispari);
				}
				for (short m_i = 0 ; m_i<dimensione_m; ++m_i) {
					parita_della_base_totale.push_back(parita_della_base_iperangolare.back());
				}
			}
		}
		dimensione_base_pari = dimensione_base_iperangolare_pari * dimensione_m;
		dimensione_base_dispari = dimensione_base_iperangolare_dispari * dimensione_m;
	
	}

	Accoppiamenti_Parita(const Accoppiamenti_Iperangolari_K_max & KK, 
						 const Accoppiamenti_Spinoriali_S & SS, 
						 const Accoppiamenti_Spinoriali_S_max & TT, 
						 const unsigned short & dimensione_m){
		dimensione_base_pari = 0;
		dimensione_base_iperangolare_pari = 0;
		dimensione_base_iperangolare_dispari = 0;
		dimensione_base_iperangolare = KK.D()*SS.D()*TT.D();
		dimensione_base_totale = dimensione_base_iperangolare*dimensione_m;
		unsigned short n_cluster = KK.P_l()[0].size();
		for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
			for (short s_i = 0 ; s_i<SS.D(); ++s_i) {
				for (short t_i = 0 ; t_i<TT.D(); ++t_i) {
					std::vector<unsigned short> p_temp;
					for (short i = 0 ; i<n_cluster; ++i) {
						p_temp.push_back((SS.P_S()[s_i][i] + KK.P_l()[k_i][i] + TT.P_S()[t_i][i])%2);
					}
					unsigned short pp = linear_sum(p_temp);
					if(pp == 0){//tutti pari
						parita_della_base_iperangolare.push_back(0);
						dimensione_base_iperangolare_pari++;
						indice_p.push_back(dimensione_base_iperangolare_pari);
						indice_d.push_back(0);
					}
					else if (pp == n_cluster){//tutti dispari
						parita_della_base_iperangolare.push_back(1);
						dimensione_base_iperangolare_dispari++;
						indice_p.push_back(0);
						indice_d.push_back(dimensione_base_iperangolare_dispari);
					}
					else{//misto
						parita_della_base_iperangolare.push_back(2);
						indice_p.push_back(0);
						indice_d.push_back(0);
					}
					for (short m_i = 0 ; m_i<dimensione_m; ++m_i) {
						parita_della_base_totale.push_back(parita_della_base_iperangolare.back());
					}
				}
			}
		}
		dimensione_base_pari = dimensione_base_iperangolare_pari * dimensione_m;
		dimensione_base_dispari = dimensione_base_iperangolare_dispari * dimensione_m;

	}
	Accoppiamenti_Parita(const Accoppiamenti_Iperangolari_K_max & KK, 
						 const Accoppiamenti_Spinoriali_S_max & SS, 
						 const Accoppiamenti_Spinoriali_S_max & TT, 
						 const unsigned short & dimensione_m){
		dimensione_base_pari = 0;
		dimensione_base_iperangolare_pari = 0;
		dimensione_base_iperangolare_dispari = 0;
		dimensione_base_iperangolare = KK.D()*SS.D()*TT.D();
		dimensione_base_totale = dimensione_base_iperangolare*dimensione_m;
		unsigned short n_cluster = KK.P_l()[0].size();
		for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
			for (short s_i = 0 ; s_i<SS.D(); ++s_i) {
				for (short t_i = 0 ; t_i<TT.D(); ++t_i) {
					std::vector<unsigned short> p_temp;
					for (short i = 0 ; i<n_cluster; ++i) {
						p_temp.push_back((SS.P_S()[s_i][i] + KK.P_l()[k_i][i] + TT.P_S()[t_i][i])%2);
					}
					unsigned short pp = linear_sum(p_temp);
					if(pp == 0){//tutti pari
						parita_della_base_iperangolare.push_back(0);
						dimensione_base_iperangolare_pari++;
						indice_p.push_back(dimensione_base_iperangolare_pari);
						indice_d.push_back(0);
					}
					else if (pp == n_cluster){//tutti dispari
						parita_della_base_iperangolare.push_back(1);
						dimensione_base_iperangolare_dispari++;
						indice_p.push_back(0);
						indice_d.push_back(dimensione_base_iperangolare_dispari);
					}
					else{//misto
						parita_della_base_iperangolare.push_back(2);
						indice_p.push_back(0);
						indice_d.push_back(0);
					}
					for (short m_i = 0 ; m_i<dimensione_m; ++m_i) {
						parita_della_base_totale.push_back(parita_della_base_iperangolare.back());
					}
				}
			}
		}
		dimensione_base_pari = dimensione_base_iperangolare_pari * dimensione_m;
		dimensione_base_dispari = dimensione_base_iperangolare_dispari * dimensione_m;
		
	}
	
	~Accoppiamenti_Parita () {
	}

	inline unsigned long D() const {return dimensione_base_totale;}	
	inline unsigned long D_p() const {return dimensione_base_pari;}	
	inline unsigned long D_d() const {return dimensione_base_dispari;}	

	inline unsigned long D(const unsigned short & PP) const {
		if(PP == 0) return dimensione_base_pari;
		else if(PP == 1) return dimensione_base_dispari;
		else return dimensione_base_totale;
	}	

	inline unsigned long D_h() const {return dimensione_base_iperangolare;}	
	inline unsigned long D_h_p() const {return dimensione_base_iperangolare_pari;}	
	inline unsigned long D_h_d() const {return dimensione_base_iperangolare_dispari;}
	
	inline unsigned long D_h(const unsigned short & PP) const {
		if(PP == 0) return dimensione_base_iperangolare_pari;
		else if(PP == 1) return dimensione_base_iperangolare_dispari;
		else return dimensione_base_iperangolare;
	}	
	
	inline ublas::symmetric_matrix<tipo_matrici> M_p(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_p_t(dimensione_base_pari,dimensione_base_pari);
		unsigned long riga1 = 0;		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_totale[riga] == 0){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_totale[colonna] == 0){
						M_p_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_p_t;
	}
	inline ublas::symmetric_matrix<tipo_matrici> M_d(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_d_t(dimensione_base_dispari,dimensione_base_dispari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_totale[riga] == 1){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_totale[colonna] == 1){
						M_d_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_d_t;
	}

	inline ublas::symmetric_matrix<tipo_matrici> M(const ublas::symmetric_matrix<tipo_matrici> & Matrice, const unsigned short & PP) const {
		if(PP == 0) return this->M_p(Matrice);
		else if(PP == 1) return this->M_d(Matrice);
		else return Matrice;
	}

	inline ublas::symmetric_matrix<tipo_matrici> M_h_p(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_p_t(dimensione_base_iperangolare_pari,dimensione_base_iperangolare_pari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_iperangolare[riga] == 0){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_iperangolare[colonna] == 0){
						M_p_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_p_t;
	}
	inline ublas::symmetric_matrix<tipo_matrici> M_h_d(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_d_t(dimensione_base_iperangolare_dispari,dimensione_base_iperangolare_dispari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_iperangolare[riga] == 1){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_iperangolare[colonna] == 1){
						M_d_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_d_t;
	}
	
	inline ublas::symmetric_matrix<tipo_matrici> M_h(const ublas::symmetric_matrix<tipo_matrici> & Matrice, const unsigned short & PP){
		if(PP == 0) return this->M_h_p(Matrice);
		else if(PP == 1) return this->M_h_d(Matrice);
		else return Matrice;
	}
	
	std::vector <V_MATRICI_SPARSE> M_p(const std::vector <V_MATRICI_SPARSE> & Matrice) const {
		std::vector <V_MATRICI_SPARSE> M_p_t;
		
		unsigned int riga,colonna,dimensione = Matrice.front().size();
		for (unsigned long i = 0 ; i < Matrice.size() ; ++i ){
			if(Matrice[i].front().nnz()){
				MATRICI_SPARSE M_temp(dimensione_base_iperangolare_pari,dimensione_base_iperangolare_pari);
				if (dimensione==2 && Matrice[i].back().nnz()) {
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2){
							riga = i2.index1();
							colonna = i2.index2();
							 if (indice_p[riga] && indice_p[colonna]) {
								tipo_matrici val = (*i2);
								 if (std::abs(val) > err_matrici) {
										M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
									}
								}
							}								
				}
				else if (dimensione==3 && Matrice[i][1].nnz() && Matrice[i].back().nnz()) {
										
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice[i][1].begin1(); j1 !=Matrice[i][1].end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice[i][1].size1()*(i2.index1())+j2.index1();
									colonna = Matrice[i][1].size2()*(i2.index2())+j2.index2();
									if (indice_p[riga] && indice_p[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
										}
									}
								}									
				}
				else if (dimensione==4 && Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice1.size1()*(i2.index1())+j2.index1();
									colonna = Matrice1.size2()*(i2.index2())+j2.index2();
									if (indice_p[riga] && indice_p[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
										}
									}
								}									
				}
				else if(Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					bool verita = true;
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					for (short j = 3 ; verita *(j < dimensione - 1); ++j) {
						verita *= Matrice[i][j].nnz();
						if (verita) Matrice1=composizione_matrici(Matrice1,Matrice[i][j]);
					}
					
					if (verita) {
						for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
							for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
								for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
									for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
										riga = Matrice1.size1()*(i2.index1())+j2.index1();
										colonna = Matrice1.size2()*(i2.index2())+j2.index2();
										if (indice_p[riga] && indice_p[colonna]) {
											tipo_matrici val = ((*i2)*(*j2));
											if (std::abs(val) > err_matrici) {
												M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
											}
										}
									}								
					}
						
				}
				if(M_temp.nnz()){
					V_MATRICI_SPARSE Matrice_out;
					Matrice_out.push_back(M_temp);
					Matrice_out.push_back(Matrice[i].back());
					M_p_t.push_back(Matrice_out);
				}
				else
					std::cerr << "matrice droppata" << std::endl;
		}
	}
			
		return M_p_t;
	}
	std::vector <V_MATRICI_SPARSE> M_d(const std::vector <V_MATRICI_SPARSE> & Matrice) const {
		std::vector <V_MATRICI_SPARSE> M_p_t;
		unsigned int riga,colonna,dimensione = Matrice.front().size();
		for (long i = 0 ; i < Matrice.size() ; ++i ){
			if(Matrice[i].front().nnz()){
				MATRICI_SPARSE M_temp(dimensione_base_iperangolare_dispari,dimensione_base_iperangolare_dispari);
				
				if (dimensione==2 && Matrice[i].back().nnz()) {
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2){
							riga = i2.index1();
							colonna = i2.index2();
							if (indice_d[riga] && indice_d[colonna]) {
								tipo_matrici val = (*i2);
								if (std::abs(val) > err_matrici) {
									M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
								}
							}
						}								
				}
				else if (dimensione==3 && Matrice[i][1].nnz() && Matrice[i].back().nnz()) {
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice[i][1].begin1(); j1 !=Matrice[i][1].end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice[i][1].size1()*(i2.index1())+j2.index1();
									colonna = Matrice[i][1].size2()*(i2.index2())+j2.index2();
									if (indice_d[riga] && indice_d[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
										}
									}
								}									
				}
				else if (dimensione==4 && Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice1.size1()*(i2.index1())+j2.index1();
									colonna = Matrice1.size2()*(i2.index2())+j2.index2();
									if (indice_d[riga] && indice_d[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
										}
									}
								}									
				}
				else if (Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()) {
					bool verita = true;
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					for (short j = 3 ; verita *(j < dimensione - 1); ++j) {
						verita *= Matrice[i][j].nnz();
						if (verita) Matrice1=composizione_matrici(Matrice1,Matrice[i][j]);
					}
					
					if (verita) {
						
						
						for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
							for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
								for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
									for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
										riga = Matrice1.size1()*(i2.index1())+j2.index1();
										colonna = Matrice1.size2()*(i2.index2())+j2.index2();
										if (indice_d[riga] && indice_d[colonna]) {
											tipo_matrici val = ((*i2)*(*j2));
											if (std::abs(val) > err_matrici) {
												M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
											}
										}
									}
					}
				}
				
				if(M_temp.nnz()){
					V_MATRICI_SPARSE Matrice_out;
					Matrice_out.push_back(M_temp);
					Matrice_out.push_back(Matrice[i].back());
					M_p_t.push_back(Matrice_out);
				}
				else
					std::cerr << "matrice droppata" << std::endl;
			}
		}
		return M_p_t;
	}

	std::vector <V_MATRICI_SPARSE> M(const std::vector <V_MATRICI_SPARSE> & Matrice, const unsigned short & PP) const {
		if(PP == 0) return this->M_p(Matrice);
		else if(PP == 1) return this->M_d(Matrice);
		else return Matrice;
	}
	
	inline VETTORI V_p(const VETTORI & Vettore_pari) const {
		VETTORI V_t(dimensione_base_totale);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 0){
				V_t(indice) = Vettore_pari(indice1);
				++indice1;
			}
			else{
				V_t(indice) = 0.;				
			}
		}
		return V_t;
	}
	inline VETTORI V_d(const VETTORI & Vettore_dispari) const {
		VETTORI V_t(dimensione_base_totale);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 1){
				V_t(indice) = Vettore_dispari(indice1);
				++indice1;
			}
			else{
				V_t(indice) = 0.;				
			}
		}
		return V_t;
	}
	inline VETTORI V_sel_p(const VETTORI & Vettore_tot) const {
		VETTORI V_pari(dimensione_base_pari);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 0){
				V_pari(indice1) = Vettore_tot(indice);
				++indice1;
			}
		}
		return V_pari;
	}
	inline VETTORI V_sel_d(const VETTORI & Vettore_tot) const {
		VETTORI V_dispari(dimensione_base_dispari);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 1){
				V_dispari(indice1) = Vettore_tot(indice);
				++indice1;
			}
		}
		return V_dispari;
	}
	
	inline VETTORI V(const VETTORI & Vettore, const unsigned short & PP){
		if(PP == 0) return this->V_p(Vettore);
		else if(PP == 1) return this->V_d(Vettore);
		else return Vettore;
	}
	
};
//TODO: eventualmente le parità miste
class Accoppiamenti_Parita_mista {
	
	typedef MATRICI_SPARSE::const_iterator1 i1_t;
	typedef MATRICI_SPARSE::const_iterator2 i2_t;
	unsigned long dimensione_base_totale;
	unsigned long dimensione_base_pari;
	unsigned long dimensione_base_dispari;
	unsigned long dimensione_base_iperangolare;
	unsigned long dimensione_base_iperangolare_pari;
	unsigned long dimensione_base_iperangolare_dispari;
	std::vector < unsigned short > parita_della_base_totale;							
	std::vector < unsigned short > parita_della_base_iperangolare;							
	std::vector < unsigned long > indice_p;							
	std::vector < unsigned long > indice_d;							
	
public:
	
	Accoppiamenti_Parita_mista(const Accoppiamenti_Iperangolari_K_max & KK, 
							   const unsigned short & dimensione_m){
		dimensione_base_pari = 0;
		dimensione_base_iperangolare_pari = 0;
		dimensione_base_iperangolare_dispari = 0;
		dimensione_base_iperangolare = KK.D();
		dimensione_base_totale = dimensione_base_iperangolare*dimensione_m;
		unsigned short n_cluster = 2;
		for (unsigned long k_i = 0 ; k_i<KK.D(); ++k_i) {
			if(KK.P_l(k_i)[0] == 0 && KK.P_l(k_i)[1] == 1 ){//tutti pari
				parita_della_base_iperangolare.push_back(0);
				dimensione_base_iperangolare_pari++;
				indice_p.push_back(dimensione_base_iperangolare_pari);
				indice_d.push_back(0);
			}
			else if (KK.P_l(k_i)[1] == 0 && KK.P_l(k_i)[0] == 1 ){//tutti dispari
				parita_della_base_iperangolare.push_back(1);
				dimensione_base_iperangolare_dispari++;
				indice_p.push_back(0);
				indice_d.push_back(dimensione_base_iperangolare_dispari);
			}
			else{//non misto
				parita_della_base_iperangolare.push_back(2);
				indice_p.push_back(0);
				indice_d.push_back(0);
			}
			for (short m_i = 0 ; m_i<dimensione_m; ++m_i) {
				parita_della_base_totale.push_back(parita_della_base_iperangolare.back());
			}
		}
		dimensione_base_pari = dimensione_base_iperangolare_pari * dimensione_m;
		dimensione_base_dispari = dimensione_base_iperangolare_dispari * dimensione_m;
	}
		
	~Accoppiamenti_Parita_mista () {
	}
	
	inline unsigned long D() const {return dimensione_base_totale;}	
	inline unsigned long D_p() const {return dimensione_base_pari;}	
	inline unsigned long D_d() const {return dimensione_base_dispari;}	
	
	inline unsigned long D(const unsigned short & PP) const {
		if(PP == 0) return dimensione_base_pari;
		else if(PP == 1) return dimensione_base_dispari;
		else return dimensione_base_totale;
	}	
	
	inline unsigned long D_h() const {return dimensione_base_iperangolare;}	
	inline unsigned long D_h_p() const {return dimensione_base_iperangolare_pari;}	
	inline unsigned long D_h_d() const {return dimensione_base_iperangolare_dispari;}
	
	inline unsigned long D_h(const unsigned short & PP) const {
		if(PP == 0) return dimensione_base_iperangolare_pari;
		else if(PP == 1) return dimensione_base_iperangolare_dispari;
		else return dimensione_base_iperangolare;
	}	
	
	inline ublas::symmetric_matrix<tipo_matrici> M_p(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_p_t(dimensione_base_pari,dimensione_base_pari);
		unsigned long riga1 = 0;		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_totale[riga] == 0){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_totale[colonna] == 0){
						M_p_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_p_t;
	}
	inline ublas::symmetric_matrix<tipo_matrici> M_d(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_d_t(dimensione_base_dispari,dimensione_base_dispari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_totale[riga] == 1){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_totale[colonna] == 1){
						M_d_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_d_t;
	}
	
	inline ublas::symmetric_matrix<tipo_matrici> M(const ublas::symmetric_matrix<tipo_matrici> & Matrice, const unsigned short & PP) const {
		if(PP == 0) return this->M_p(Matrice);
		else if(PP == 1) return this->M_d(Matrice);
		else return Matrice;
	}
	
	inline ublas::symmetric_matrix<tipo_matrici> M_h_p(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_p_t(dimensione_base_iperangolare_pari,dimensione_base_iperangolare_pari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_iperangolare[riga] == 0){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_iperangolare[colonna] == 0){
						M_p_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_p_t;
	}
	inline ublas::symmetric_matrix<tipo_matrici> M_h_d(const ublas::symmetric_matrix<tipo_matrici> & Matrice) const {
		ublas::symmetric_matrix<tipo_matrici> M_d_t(dimensione_base_iperangolare_dispari,dimensione_base_iperangolare_dispari);
		unsigned long riga1 = 0;
		
		for (long riga = 0; riga < Matrice.size1(); ++riga) {
			if (parita_della_base_iperangolare[riga] == 1){
				unsigned long colonna1 = 0;
				for (long colonna = 0; colonna <= riga; ++colonna) {
					if (parita_della_base_iperangolare[colonna] == 1){
						M_d_t(riga1,colonna1) = Matrice(riga,colonna);
						++colonna1;
					}
				}
				++riga1;
			}
		}
		return M_d_t;
	}
	
	inline ublas::symmetric_matrix<tipo_matrici> M_h(const ublas::symmetric_matrix<tipo_matrici> & Matrice, const unsigned short & PP){
		if(PP == 0) return this->M_h_p(Matrice);
		else if(PP == 1) return this->M_h_d(Matrice);
		else return Matrice;
	}
	
	std::vector <V_MATRICI_SPARSE> M_p(const std::vector <V_MATRICI_SPARSE> & Matrice) const {
		std::vector <V_MATRICI_SPARSE> M_p_t;
		
		unsigned int riga,colonna,dimensione = Matrice.front().size();
		for (unsigned long i = 0 ; i < Matrice.size() ; ++i ){
			if(Matrice[i].front().nnz()){
				MATRICI_SPARSE M_temp(dimensione_base_iperangolare_pari,dimensione_base_iperangolare_pari);
				if (dimensione==2 && Matrice[i].back().nnz()) {
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2){
							riga = i2.index1();
							colonna = i2.index2();
							if (indice_p[riga] && indice_p[colonna]) {
								tipo_matrici val = (*i2);
								if (std::abs(val) > err_matrici) {
									M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
								}
							}
						}								
				}
				else if (dimensione==3 && Matrice[i][1].nnz() && Matrice[i].back().nnz()) {
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice[i][1].begin1(); j1 !=Matrice[i][1].end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice[i][1].size1()*(i2.index1())+j2.index1();
									colonna = Matrice[i][1].size2()*(i2.index2())+j2.index2();
									if (indice_p[riga] && indice_p[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
										}
									}
								}									
				}
				else if (dimensione==4 && Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice1.size1()*(i2.index1())+j2.index1();
									colonna = Matrice1.size2()*(i2.index2())+j2.index2();
									if (indice_p[riga] && indice_p[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
										}
									}
								}									
				}
				else if(Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					bool verita = true;
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					for (short j = 3 ; verita *(j < dimensione - 1); ++j) {
						verita *= Matrice[i][j].nnz();
						if (verita) Matrice1=composizione_matrici(Matrice1,Matrice[i][j]);
					}
					
					if (verita) {
						for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
							for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
								for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
									for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
										riga = Matrice1.size1()*(i2.index1())+j2.index1();
										colonna = Matrice1.size2()*(i2.index2())+j2.index2();
										if (indice_p[riga] && indice_p[colonna]) {
											tipo_matrici val = ((*i2)*(*j2));
											if (std::abs(val) > err_matrici) {
												M_temp(indice_p[riga]-1,indice_p[colonna]-1) = val;
											}
										}
									}								
					}
					
				}
				if(M_temp.nnz()){
					V_MATRICI_SPARSE Matrice_out;
					Matrice_out.push_back(M_temp);
					Matrice_out.push_back(Matrice[i].back());
					M_p_t.push_back(Matrice_out);
				}
				else
					std::cerr << "matrice droppata" << std::endl;
			}
		}
		
		return M_p_t;
	}
	std::vector <V_MATRICI_SPARSE> M_d(const std::vector <V_MATRICI_SPARSE> & Matrice) const {
		std::vector <V_MATRICI_SPARSE> M_p_t;
		unsigned int riga,colonna,dimensione = Matrice.front().size();
		for (long i = 0 ; i < Matrice.size() ; ++i ){
			if(Matrice[i].front().nnz()){
				MATRICI_SPARSE M_temp(dimensione_base_iperangolare_dispari,dimensione_base_iperangolare_dispari);
				
				if (dimensione==2 && Matrice[i].back().nnz()) {
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2){
							riga = i2.index1();
							colonna = i2.index2();
							if (indice_d[riga] && indice_d[colonna]) {
								tipo_matrici val = (*i2);
								if (std::abs(val) > err_matrici) {
									M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
								}
							}
						}								
				}
				else if (dimensione==3 && Matrice[i][1].nnz() && Matrice[i].back().nnz()) {
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice[i][1].begin1(); j1 !=Matrice[i][1].end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice[i][1].size1()*(i2.index1())+j2.index1();
									colonna = Matrice[i][1].size2()*(i2.index2())+j2.index2();
									if (indice_d[riga] && indice_d[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
										}
									}
								}									
				}
				else if (dimensione==4 && Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()){
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					
					for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
						for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
							for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
								for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
									riga = Matrice1.size1()*(i2.index1())+j2.index1();
									colonna = Matrice1.size2()*(i2.index2())+j2.index2();
									if (indice_d[riga] && indice_d[colonna]) {
										tipo_matrici val = ((*i2)*(*j2));
										if (std::abs(val) > err_matrici) {
											M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
										}
									}
								}									
				}
				else if (Matrice[i][1].nnz() && Matrice[i][2].nnz() && Matrice[i].back().nnz()) {
					bool verita = true;
					MATRICI_SPARSE Matrice1(composizione_matrici(Matrice[i][1],Matrice[i][2]));
					for (short j = 3 ; verita *(j < dimensione - 1); ++j) {
						verita *= Matrice[i][j].nnz();
						if (verita) Matrice1=composizione_matrici(Matrice1,Matrice[i][j]);
					}
					
					if (verita) {
						
						
						for (i1_t i1 = Matrice[i].front().begin1(); i1 !=Matrice[i].front().end1(); ++i1)
							for (i1_t j1 = Matrice1.begin1(); j1 !=Matrice1.end1(); ++j1)
								for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
									for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
										riga = Matrice1.size1()*(i2.index1())+j2.index1();
										colonna = Matrice1.size2()*(i2.index2())+j2.index2();
										if (indice_d[riga] && indice_d[colonna]) {
											tipo_matrici val = ((*i2)*(*j2));
											if (std::abs(val) > err_matrici) {
												M_temp(indice_d[riga]-1,indice_d[colonna]-1) = val;
											}
										}
									}
					}
				}
				
				if(M_temp.nnz()){
					V_MATRICI_SPARSE Matrice_out;
					Matrice_out.push_back(M_temp);
					Matrice_out.push_back(Matrice[i].back());
					M_p_t.push_back(Matrice_out);
				}
				else
					std::cerr << "matrice droppata" << std::endl;
			}
		}
		return M_p_t;
	}
	
	std::vector <V_MATRICI_SPARSE> M(const std::vector <V_MATRICI_SPARSE> & Matrice, const unsigned short & PP) const {
		if(PP == 0) return this->M_p(Matrice);
		else if(PP == 1) return this->M_d(Matrice);
		else return Matrice;
	}
	
	inline VETTORI V_p(const VETTORI & Vettore_pari) const {
		VETTORI V_t(dimensione_base_totale);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 0){
				V_t(indice) = Vettore_pari(indice1);
				++indice1;
			}
			else{
				V_t(indice) = 0.;				
			}
		}
		return V_t;
	}
	inline VETTORI V_d(const VETTORI & Vettore_dispari) const {
		VETTORI V_t(dimensione_base_totale);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 1){
				V_t(indice) = Vettore_dispari(indice1);
				++indice1;
			}
			else{
				V_t(indice) = 0.;				
			}
		}
		return V_t;
	}
	inline VETTORI V_sel_p(const VETTORI & Vettore_tot) const {
		VETTORI V_pari(dimensione_base_pari);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 0){
				V_pari(indice1) = Vettore_tot(indice);
				++indice1;
			}
		}
		return V_pari;
	}
	inline VETTORI V_sel_d(const VETTORI & Vettore_tot) const {
		VETTORI V_dispari(dimensione_base_dispari);
		unsigned long indice1 = 0;
		for (long indice = 0; indice < dimensione_base_totale; ++indice) {
			if (parita_della_base_totale[indice] == 1){
				V_dispari(indice1) = Vettore_tot(indice);
				++indice1;
			}
		}
		return V_dispari;
	}
	
	inline VETTORI V(const VETTORI & Vettore, const unsigned short & PP){
		if(PP == 0) return this->V_p(Vettore);
		else if(PP == 1) return this->V_d(Vettore);
		else return Vettore;
	}
	
};
#endif
