/*
 *  inversione_matrice_ublas.h
 *  Tesi
 *
 *  Created by Fredrik Orderud http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
 *
 */

#ifndef INVERT_MATRIX_HPP
#define INVERT_MATRIX_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp> 
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/operation_sparse.hpp> 
#include <boost/numeric/ublas/lu.hpp>
#include <assert.h>
#include "tipi_definiti.h"
#include "combinazioni.h"
#include <omp.h> //calcolo parallelo OMP
#include <boost/timer.hpp>
#include <mkl_service.h>
#include "mkl_prod_marco.h"
#include "err.h"




namespace boost { namespace numeric { namespace ublas {
	
	inline VETTORI sqrt_vec(const VETTORI & V){
		VETTORI V_temp(V.size());
		for (short i = 0; i<V.size(); ++i) {
			V_temp(i) = sqrt(V(i));
		}
		return V_temp;
	}

	inline long double quad_sum_vec(const VETTORI & V){
		long double res = 0.;
#pragma omp parallel for reduction(+:res) schedule( guided )
		for (short i = 0; i<V.size(); ++i) 
			res += pow(V(i),2);
		
		return res;
	}

	inline long double quad_sum_mat(const MATRICI & M){
		long double res = 0.;
#pragma omp parallel for reduction(+:res) schedule( guided )
		for (short i = 0; i<M.size1(); ++i) 
			for (short j = 0; j<M.size2(); ++j) 
			res += pow(M(i,j),2);

		return res;
	}
	
	inline MATRICI_SIMMETRICHE_UP sqrt_mat_sim_up(const MATRICI_SIMMETRICHE_UP & M){
		MATRICI_SIMMETRICHE_UP M_temp(M.size1(),M.size2());
		for (short i = 0; i<M.size1(); ++i) {
			for (short j = i; j < M.size2(); ++j) {
				M_temp(i,j) = sqrt(M(i,j));
			}
		}
		return M_temp;
	}
	
	
	inline MATRICI_SPARSE composizione_matrici(const MATRICI_SPARSE & M1,
											   const MATRICI_SPARSE & M2){
		//boost::timer t;
		assert(M1.size1() == M1.size2());
		assert(M2.size1() == M2.size2());
		unsigned long D_tot = M1.size1()*M2.size1(),riga,colonna;
		tipo_matrici val;
		MATRICI_SPARSE M(D_tot,D_tot);
		typedef MATRICI_SPARSE::const_iterator1 i1_t;
		typedef MATRICI_SPARSE::const_iterator2 i2_t;
		
		for (i1_t i1 = M1.begin1(); i1 !=M1.end1(); ++i1)
			for (i1_t j1 = M2.begin1(); j1 !=M2.end1(); ++j1)
				for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
					for (i2_t j2 = j1.begin(); j2 != j1.end(); ++j2){
						riga = M2.size1()*(i2.index1())+j2.index1();
						colonna = M2.size2()*(i2.index2())+j2.index2();
						val = ((*i2)*(*j2));
						if (std::abs(val) > errore_matrici) {
							M.push_back(riga,colonna,val);
						}
					}
		return M;
	}

	
	inline symmetric_matrix<tipo_matrici> composizione_matrici_sim(const V_MATRICI_SPARSE & M){
		BOOST_UBLAS_CHECK (M[0].size1() == M[0].size2(), bad_size ()); //quadrate
		MATRICI_SPARSE M_o(M.front());
		for (short i = 1; i<M.size(); ++i) {
			BOOST_UBLAS_CHECK (M[i].size1() == M[i].size2(), bad_size ()); //quadrate
			M_o = composizione_matrici(M_o,M[i]);
		}
		return M_o;
	} 
	
	inline MATRICI_SPARSE composizione_matrici(const V_MATRICI_SPARSE & M){
		BOOST_UBLAS_CHECK (M[0].size1() == M[0].size2(), bad_size ()); //quadrate
		MATRICI_SPARSE M_o(M.front());
		for (short i = 1; i<M.size(); ++i) {
			BOOST_UBLAS_CHECK (M[i].size1() == M[i].size2(), bad_size ()); //quadrate
			M_o = composizione_matrici(M_o,M[i]);
		}
		return M_o;
	} 
	
	symmetric_matrix<tipo_matrici> composizione_matrici(const V_V_MATRICI_SPARSE & M){
		symmetric_matrix<tipo_matrici> M_o(composizione_matrici(M[0]));
		for (short i = 1; i<M.size(); ++i) {
			if (M[i].front().nnz()) {
				BOOST_UBLAS_CHECK (M[i].size1() == M[i].size2(), bad_size ()); //quadrate
				symmetric_matrix<tipo_matrici> M_temp(composizione_matrici(M[i]));
				M_o.plus_assign(M_temp);				
			}
		}
		return M_o;
	} 
		
	
	
	MATRICI_SPARSE accresci( const MATRICI_SPARSE & m, const unsigned short & accrescimento) { 
	
		MATRICI_SPARSE A(m.size1() + accrescimento, m.size2() + accrescimento);

		typedef MATRICI_SPARSE::const_iterator1 i1_t;
		typedef MATRICI_SPARSE::const_iterator2 i2_t;
	
		for (i1_t i1 = m.begin1(); i1 !=m.end1(); ++i1)
			for (i2_t i2 = i1.begin(); i2 != i1.end(); ++i2)
				A.push_back(i2.index1(),i2.index2(),*i2);
			
		return A;
	}

	
	
	/** General matrix determinant. 
	 * It uses lu_factorize in uBLAS. 
	 */ 
	template<class T> 
	double lu_det(const matrix<T> & m) { 
        matrix<T> mLu(m); 
        ublas::permutation_matrix<std::size_t> pivots(m.size1()); 
        lu_factorize(mLu, pivots); 
        double det = 1.0; 
        for (std::size_t i=0; i < pivots.size(); ++i) { 
			if (pivots(i) != i) 
				det *= -1.0; 
			det *= mLu(i,i); 
        } 
        return det; 
	} 
	
	
	/* Matrix inversion routine.
	 *
	 *  Created by Fredrik Orderud http://www.crystalclearsoftware.com/cgi-bin/boost_wiki/wiki.pl?LU_Matrix_Inversion
	 *	Uses lu_factorize and lu_substitute in uBLAS to invert a matrix 
	 */
	
	
	template<class T>
	bool inverti_matrice (const matrix<T>& input, matrix<T>& inverse) {
//		typedef permutation_matrix<std::size_t> pmatrix;
		// create a working copy of the input
		matrix<T> A(input);
		// create a permutation matrix for the LU-factorization
		permutation_matrix<std::size_t> pm(A.size1());
		// perform LU-factorization
		int res = lu_factorize(A,pm);
		if( res != 0 ) 
			return false;
		// create identity matrix of "inverse"
		inverse.assign(ublas::identity_matrix<T>(A.size1()));
		// backsubstitute to get the inverse
		lu_substitute(A, pm, inverse);
		return true;
	}


	
	//Matrice inversa
	template<class T>
	matrix<T> inversa (const matrix<T>& input) {
		matrix<T> inverse(input);
		bool z = inverti_matrice(input, inverse);
		assert(z == true);	
		return inverse;
	}


/*
 
 
	//vettore (m*n) da matrice densa (m x n)
	template<class T>
	inline vector<T> m2v(const matrix<T> & V){	
		vector<T> v (V.size1()*V.size2());
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j) 
				v(i+V.size1()*j)=V(i,j);
		return v;
	}

	//matrice densa (m x n) da vettore (m*n)
	template<class T>
	inline matrix<T> v2m(const vector<T> & v,
						 const unsigned short & dim_base1){	
		assert(v.size()%dim_base1 == 0);
		unsigned short dim_base2 = v.size()/dim_base1;
		matrix<T> V (dim_base2,dim_base1);
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j) 
				V(i,j) = v(i+V.size1()*j);
		return V;
	}

	
	inline MATRICI v2m_2(const VETTORI & v,
						 const unsigned short & dim_base1){	
		assert(v.size()%dim_base1 == 0);
		unsigned short dim_base2 = v.size()/dim_base1;
		MATRICI V (dim_base2,dim_base1);
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j) 
				V(i,j) = v(i+V.size1()*j);
		return V;
	}
	
	//da un vettore genero un vettore di matrici che segue la mia formula
	
	inline std::vector<MATRICI> v2m_3(const VETTORI & v,
									  const vettore_ushort & dim_base) {
		assert(dim_base.size()==3);
		MATRICI M1(v2m_2(v,dim_base[0]));
		std::vector<MATRICI> M2;
		for (short i = 0; i<dim_base[0]; ++i) {
			M2.push_back(v2m_2((VETTORI)column(M1, i),dim_base[1]));
		}
		return M2;
	}
	
	inline std::vector<V_MATRICI> v2m_4(const VETTORI & v,
										const vettore_ushort & dim_base) {
		assert(dim_base.size()==4);
		MATRICI M1(v2m_2(v,dim_base[0]));

		//trasferisco da matrice a vettore di vettori colonna D=2: l'indice di vettore è nello spazio 0
		std::vector<V_MATRICI> M3;
		for (short i = 0; i<dim_base[0]; ++i) {
			MATRICI M2(v2m_2((VETTORI)column(M1, i),dim_base[1]));
			std::vector<MATRICI> M3_t;
			for (short j = 0; j<dim_base[1]; ++j) {
				M3_t.push_back(v2m_2((VETTORI)column(M2, j),dim_base[2]));
			}
			M3.push_back(M3_t);
		}
		return M3;
	}
*/
	
	
	
	//vettore (m*n) da matrice densa (m x n)
	template<class T>
	inline vector<T> m2v(const matrix<T> & V){	
		vector<T> v (V.size1()*V.size2());
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j){
				if (std::abs(V(i,j)) > errore_matrici) {
					v(j+V.size2()*i)=V(i,j);
				}
				else {
					v(j+V.size2()*i)=0.;
				}

			}
				
		return v;
	}
	
	//matrice densa (m x n) da vettore (m*n)
	template<class T>
	inline matrix<T> v2m(const vector<T> & v,
						 const unsigned short & dim_base1){	
		assert(v.size()%dim_base1 == 0);
		unsigned short dim_base2 = v.size()/dim_base1;
		matrix<T> V (dim_base1,dim_base2);
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j) {
				if (std::abs(v(j+V.size2()*i)) > errore_matrici) {
					V(i,j) = v(j+V.size2()*i);
				}
				else {
					V(i,j) = 0.;
				}
				
			}
		return V;
	}
	
	
	inline MATRICI v2m_2(const VETTORI & v,
						 const unsigned short & dim_base1){	
		assert(v.size()%dim_base1 == 0);
		unsigned short dim_base2 = v.size()/dim_base1;
		MATRICI V (dim_base1,dim_base2);
		for (short i = 0; i < V.size1(); ++i) 
			for (short j = 0; j < V.size2(); ++j) {
				if (std::abs(v(j+V.size2()*i)) > errore_matrici) {
					V(i,j) = v(j+V.size2()*i);
				}
				else {
					V(i,j) = 0.;
				}
				
			}
		return V;
	}
	
	//da un vettore genero un vettore di matrici che segue la mia formula
	
	inline std::vector<MATRICI> v2m_3(const VETTORI & v,
									  const vettore_ushort & dim_base) {
		assert(dim_base.size()==3);
		MATRICI M1(v2m_2(v,dim_base[0]));
		std::vector<MATRICI> M2;
		for (short i = 0; i<dim_base[0]; ++i) {
			M2.push_back(v2m_2((VETTORI)row(M1, i),dim_base[1]));
		}
		return M2;
	}
	
	inline std::vector<V_MATRICI> v2m_4(const VETTORI & v,
										const vettore_ushort & dim_base) {
		assert(dim_base.size()==4);
		MATRICI M1(v2m_2(v,dim_base[0]));
		
		//trasferisco da matrice a vettore di vettori colonna D=2: l'indice di vettore è nello spazio 0
		std::vector<V_MATRICI> M3;
		for (short i = 0; i<dim_base[0]; ++i) {
			MATRICI M2(v2m_2((VETTORI)row(M1, i),dim_base[1]));
			std::vector<MATRICI> M3_t;
			for (short j = 0; j<dim_base[1]; ++j) {
				M3_t.push_back(v2m_2((VETTORI)row(M2, j),dim_base[2]));
			}
			M3.push_back(M3_t);
		}
		return M3;
	}	

/*	inline MATRICI mult_2(const V_MATRICI_SPARSE & M,
						  const MATRICI & x) {
		//		return prod_mkl_sparse_trans(M[0],(MATRICI)trans(prod_mkl((MATRICI)M[1],x))); ma il secondo pezzo è X' * M[1] (che è simmetrica)
		return prod_mkl_sparse_trans(M[0],prod_mkl_trans1( trans ( x ),(MATRICI)M[1]));
	}
*/
/*	inline MATRICI mult_2(const V_MATRICI_SPARSE & M,
						  const MATRICI & x) {
		//		return prod_mkl( prod_mkl((MATRICI)M[1],x),(MATRICI)M[0]);	
		return prod_mkl_sparse_trans(M[0],(MATRICI)trans(prod_mkl((MATRICI)M[1],x)));
	}
*/
	inline MATRICI mult_2(const V_MATRICI_SPARSE & M,
						  const MATRICI & x) {
		//		return prod_mkl( prod_mkl((MATRICI)M[1],x),(MATRICI)M[0]);	
		return prod_mkl_sparse(M[0],(prod_mkl(x,(MATRICI)M[1])));
	}
	
	inline MATRICI mult_2(const MATRICI_SPARSE & M_t_0,
						  const MATRICI_SPARSE & M_t_1,
						  const MATRICI & x) {
		return prod_mkl( prod_mkl((MATRICI)M_t_0,x),(MATRICI)M_t_1);		
	}

	inline void mult_2(const V_MATRICI_SPARSE & M,
					   const MATRICI & x,
					   MATRICI & V) {
		//		return prod_mkl( prod_mkl((MATRICI)M[1],x),(MATRICI)M[0]);	
		prod_mkl_sparse(M[0],(prod_mkl(x,(MATRICI)M[1])),V);
	}
	
	inline MATRICI mult_3(const V_MATRICI_SPARSE & M,
						  const std::vector<MATRICI> & x) {
		
		//mi ricavo le dimensioni dello spazio vettoriale
		vettore_ushort dim_base;
		short num_spazi = M.size();
		BOOST_UBLAS_CHECK (num_spazi == 3, bad_size ()); //quadrate
		unsigned int D=1;
		for (short i = 0; i<num_spazi; ++i) {
			BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
			dim_base.push_back(M[i].size1());
			D *=dim_base.back();
		}
		BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni	
		
		MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]);
		for ( long i = 0; i < dim_base[0]; ++i) {			
			row(M2, i) = m2v(prod_mkl(prod_mkl((MATRICI)M[1],x[i]),(MATRICI)M[2]));
		}
		return prod_mkl_sparse(M[0],M2);		
	}
	
	inline MATRICI mult_3(const MATRICI_SPARSE & M_t_0,
						  const MATRICI_SPARSE & M_t_1,
						  const MATRICI_SPARSE & M_t_2,
						  const std::vector<MATRICI> & x) {
		
		//mi ricavo le dimensioni dello spazio vettoriale
		vettore_ushort dim_base;
		dim_base.push_back(M_t_0.size1());
		dim_base.push_back(M_t_1.size1());
		dim_base.push_back(M_t_2.size1());
		
		MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]);
		for (long i = 0; i < dim_base[0]; ++i) {
			row(M2, i) = m2v(prod_mkl(prod_mkl((MATRICI)M_t_1,x[i]),(MATRICI)M_t_2));
		}
		
		return prod_mkl((MATRICI)M_t_0,M2);		
	}

	inline void mult_3(const V_MATRICI_SPARSE & M,
					   const std::vector<MATRICI> & x,
					   MATRICI & V) {
		
		//mi ricavo le dimensioni dello spazio vettoriale
		vettore_ushort dim_base;
		short num_spazi = M.size();
		BOOST_UBLAS_CHECK (num_spazi == 3, bad_size ()); //quadrate
		unsigned int D=1;
		for (short i = 0; i<num_spazi; ++i) {
			BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
			dim_base.push_back(M[i].size1());
			D *=dim_base.back();
		}
		BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni	
		
		MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]);
#pragma omp parallel for schedule( guided )
		for ( long i = 0; i < dim_base[0]; ++i) {			
			row(M2, i) = m2v(prod_mkl(prod_mkl((MATRICI)M[1],x[i]),(MATRICI)M[2]));
		}
		return prod_mkl_sparse(M[0],M2,V);		
	}
	
	
	inline MATRICI mult_4(const V_MATRICI_SPARSE & M,
						  const std::vector<V_MATRICI> & x) {
		vettore_ushort dim_base;
		short num_spazi = M.size();
		BOOST_UBLAS_CHECK (num_spazi == 4, bad_size ()); //quadrate
		unsigned int D=1;
		for (short i = 0; i<num_spazi; ++i) {
			BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
			dim_base.push_back(M[i].size1());
			D *=dim_base.back();
		}
		BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni
		
		MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]*dim_base[3]);
		
		for (long i = 0; i < dim_base[0]; ++i) {
			
			MATRICI M2_b(dim_base[1],dim_base[2]*dim_base[3]);
			
		#pragma omp parallel for schedule( guided )
			for (long i1 = 0; i1 < dim_base[1]; ++i1) {
				row(M2_b, i1) = m2v(prod_mkl(prod_mkl((MATRICI)M[2],x[i][i1]),(MATRICI)M[3]));
			}
			row(M2, i) = m2v(prod_mkl((MATRICI)M[1],M2_b));
		}	
		
		return prod_mkl_sparse(M[0],M2);		
		
	}

	inline void mult_4(const V_MATRICI_SPARSE & M,
						  const std::vector<V_MATRICI> & x,
						  MATRICI & V) {
		vettore_ushort dim_base;
		short num_spazi = M.size();
		BOOST_UBLAS_CHECK (num_spazi == 4, bad_size ()); //quadrate
		unsigned int D=1;
		for (short i = 0; i<num_spazi; ++i) {
			BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
			dim_base.push_back(M[i].size1());
			D *=dim_base.back();
		}
		BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni
		
		MATRICI M2(dim_base[0],dim_base[1]*dim_base[2]*dim_base[3]);
		
		for (long i = 0; i < dim_base[0]; ++i) {
			
			MATRICI M2_b(dim_base[1],dim_base[2]*dim_base[3]);
			
#pragma omp parallel for schedule( guided )
			for (long i1 = 0; i1 < dim_base[1]; ++i1) {
				row(M2_b, i1) = m2v(prod_mkl(prod_mkl((MATRICI)M[2],x[i][i1]),(MATRICI)M[3]));
			}
			row(M2, i) = m2v(prod_mkl((MATRICI)M[1],M2_b));
		}	
		
		return prod_mkl_sparse(M[0],M2,V);		
		
	}
	
	VETTORI v_mult(const V_V_MATRICI_SPARSE & M,
				   const VETTORI & x) {
		vettore_ushort dim_base;
		for (short i = 0; i<M[0].size(); ++i) {
			dim_base.push_back(M[0][i].size1());
		}
		if (dim_base.size() == 2) {
			MATRICI V(v2m_2(x,dim_base[0]));
			//			std::cerr << " MATRICE 0" << std::endl;
			MATRICI V_o(mult_2(M[0],V));
			for (short i = 1; i < M.size(); ++i) {
				//	std::cerr << " MATRICE " << i << std::endl;
				if (M[i].front().nnz()) mult_2(M[i],V,V_o);
				//V_o.plus_assign(mult_2(M[i],V));
			}
			mkl_free_buffers();
			return m2v(V_o);			
		}
		else if (dim_base.size() == 3) {
			V_MATRICI V(v2m_3(x,dim_base));
			MATRICI V_o(mult_3(M[0],V));
			for (short i = 1; i < M.size(); ++i) {
				if (M[i].front().nnz()) mult_3(M[i],V,V_o);
				//V_o.plus_assign(mult_3(M[i],V));
			}
			mkl_free_buffers();
			return m2v(V_o);			
		}
		else if (dim_base.size() == 4) {
			V_V_MATRICI V(v2m_4(x,dim_base));
			MATRICI V_o(mult_4(M[0],V));
			for (short i = 1; i < M.size(); ++i) {
				if (M[i].front().nnz()) mult_4(M[i],V,V_o);
					//V_o.plus_assign(mult_4(M[i],V));
			}
			mkl_free_buffers();
			return m2v(V_o);			
		}
		else return x;
	}
}}}

/*
inline MATRICI mult_2(const V_MATRICI_SPARSE & M,
					  const MATRICI & x) {
	return prod_mkl(MATRICI (prod_mkl_sparse(M[1],x)),(MATRICI)M[0]);		
}

inline MATRICI mult_2(const MATRICI_SPARSE & M_t_0,
					  const MATRICI_SPARSE & M_t_1,
					  const MATRICI & x) {
	return prod_mkl(MATRICI (prod_mkl_sparse(M_t_1,x)),(MATRICI)M_t_0);		
}

inline MATRICI mult_3(const V_MATRICI_SPARSE & M,
					  const std::vector<MATRICI> & x) {
	
	//mi ricavo le dimensioni dello spazio vettoriale
	vettore_ushort dim_base;
	short num_spazi = M.size();
	BOOST_UBLAS_CHECK (num_spazi == 3, bad_size ()); //quadrate
	unsigned int D=1;
	for (short i = 0; i<num_spazi; ++i) {
		BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
		dim_base.push_back(M[i].size1());
		D *=dim_base.back();
	}
	BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni	
	
	MATRICI M2(dim_base[1]*dim_base[2],dim_base[0]);
	for ( long i = 0; i < dim_base[0]; ++i) {			
		column(M2, i) = m2v(prod_mkl(MATRICI (prod_mkl_sparse(M[2],x[i])),(MATRICI)M[1]));
	}
	return prod_mkl(M2,(MATRICI)M[0]);		
}

inline MATRICI mult_3(const MATRICI_SPARSE & M_t_0,
					  const MATRICI_SPARSE & M_t_1,
					  const MATRICI_SPARSE & M_t_2,
					  const std::vector<MATRICI> & x) {
	
	//mi ricavo le dimensioni dello spazio vettoriale
	vettore_ushort dim_base;
	dim_base.push_back(M_t_0.size1());
	dim_base.push_back(M_t_1.size1());
	dim_base.push_back(M_t_2.size1());
	
	MATRICI M2(dim_base[1]*dim_base[2],dim_base[0]);
	for (long i = 0; i < dim_base[0]; ++i) {
		column(M2, i) = m2v(prod_mkl(MATRICI (prod_mkl_sparse(M_t_2,x[i])),(MATRICI)M_t_1));
	}
	
	return prod_mkl(M2,(MATRICI)M_t_0);		
}




inline MATRICI mult_4(const V_MATRICI_SPARSE & M,
					  const std::vector<V_MATRICI> & x) {
	vettore_ushort dim_base;
	short num_spazi = M.size();
	BOOST_UBLAS_CHECK (num_spazi == 4, bad_size ()); //quadrate
	unsigned int D=1;
	for (short i = 0; i<num_spazi; ++i) {
		BOOST_UBLAS_CHECK (M[i].size1()==M[i].size2(), bad_size ()); //quadrate
		dim_base.push_back(M[i].size1());
		D *=dim_base.back();
	}
	BOOST_UBLAS_CHECK (D==x.size(), bad_size ()); //vettori di giuste dimensioni
	
	MATRICI M2(dim_base[1]*dim_base[2]*dim_base[3],dim_base[0]);
	
	for (long i = 0; i < dim_base[0]; ++i) {
		
		MATRICI M2_b(dim_base[2]*dim_base[3],dim_base[1]);
		
		for (long i1 = 0; i1 < dim_base[1]; ++i1) {
			column(M2_b, i1) = m2v(prod_mkl(MATRICI (prod_mkl_sparse(M[3],x[i][i1])),(MATRICI)M[2]));
		}
		column(M2, i) = m2v(prod_mkl(M2_b,(MATRICI)M[1]));
	}	
	
	return prod_mkl(M2,(MATRICI)M[0]);		
	
}

*/


#endif 