/*
 *  mkl_prod_marco.h
 *  Tesi
 *
 *  Created by Marco Resa on 15/03/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef _MKL_BOOST_UBLAS_MATRIX_PROD_MARCO_
#define _MKL_BOOST_UBLAS_MATRIX_PROD_MARCO_

#ifdef NDEBUG

#include <boost/version.hpp>
#if defined (BOOST_VERSION) && \
((BOOST_VERSION == 103401) \
|| (BOOST_VERSION == 103500) \
|| (BOOST_VERSION == 103600) \
|| (BOOST_VERSION == 103700) \
|| (BOOST_VERSION == 103800) \
|| (BOOST_VERSION == 103900) \
|| (BOOST_VERSION == 104000) \
|| (BOOST_VERSION == 104100) \
|| (BOOST_VERSION == 104200))

#include <boost/numeric/ublas/detail/concepts.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp> 
#include <mkl_cblas.h>
#include <mkl_spblas.h>


#ifndef MKL_BOOST_UBLAS_INLINE
#define MKL_BOOST_UBLAS_INLINE inline
#endif
namespace boost { namespace numeric { namespace ublas { 
	namespace mkl {
		
		template<class E1, class E2>
		BOOST_UBLAS_INLINE
		typename matrix_matrix_binary_traits<typename E1::value_type, E1,
		typename E2::value_type, E2>::result_type
		gemm(const matrix_expression<E1> &e1,
			 const matrix_expression<E2> &e2)
		{
			typedef typename matrix_matrix_binary_traits<typename E1::value_type, E1,
			typename E2::value_type, E2>::storage_category storage_category;
			typedef typename matrix_matrix_binary_traits<typename E1::value_type, E1,
			typename E2::value_type, E2>::orientation_category orientation_category;
			return prod (e1, e2, storage_category (), orientation_category ());
		} // For unsupported matrix types.
		template<class T>
		MKL_BOOST_UBLAS_INLINE
		void gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
				  const  CBLAS_TRANSPOSE TransB, const MKL_INT m, const MKL_INT n,
				  const MKL_INT k, const T alpha, const T *a,
				  const MKL_INT lda, const T *b, const MKL_INT ldb,
				  const T beta, T *c, const MKL_INT ldc)
		{}// To resolve externals for unsupported matrix types. Never is called.
		template<>
		MKL_BOOST_UBLAS_INLINE
		void
		gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
			 const  CBLAS_TRANSPOSE TransB, const MKL_INT m, const MKL_INT n,
			 const MKL_INT k, const double alpha, const double *a,
			 const MKL_INT lda, const double *b, const MKL_INT ldb,
			 const double beta, double *c, const MKL_INT ldc)
		{
			cblas_dgemm(Order, TransA, TransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
		}
		template<>
		MKL_BOOST_UBLAS_INLINE
		void
		gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
			 const  CBLAS_TRANSPOSE TransB, const MKL_INT m, const MKL_INT n,
			 const MKL_INT k, const float alpha, const float *a,
			 const MKL_INT lda, const float *b, const MKL_INT ldb,
			 const float beta, float *c, const MKL_INT ldc)
		{
			cblas_sgemm(Order, TransA, TransB, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);
		}
		template<>
		MKL_BOOST_UBLAS_INLINE
		void
		gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
			 const  CBLAS_TRANSPOSE TransB, const MKL_INT m, const MKL_INT n,
			 const MKL_INT k, const std::complex<double> alpha, const std::complex<double> *a,
			 const MKL_INT lda, const std::complex<double> *b, const MKL_INT ldb,
			 const std::complex<double> beta, std::complex<double> *c, const MKL_INT ldc)
		{
			cblas_zgemm(Order, TransA, TransB, m, n, k, static_cast<const void *>(&alpha),
						static_cast<const void *>(a), lda, static_cast<const void *>(b), ldb,
						static_cast<const void *>(&beta), static_cast<void *>(c), ldc);
		}
		template<>
		MKL_BOOST_UBLAS_INLINE
		void
		gemm(const  CBLAS_ORDER Order, const  CBLAS_TRANSPOSE TransA,
			 const  CBLAS_TRANSPOSE TransB, const MKL_INT m, const MKL_INT n,
			 const MKL_INT k, const std::complex<float> alpha, const std::complex<float> *a,
			 const MKL_INT lda, const std::complex<float> *b, const MKL_INT ldb,
			 const std::complex<float> beta, std::complex<float> *c, const MKL_INT ldc)
		{
			cblas_cgemm(Order, TransA, TransB, m, n, k, static_cast<const void *>(&alpha),
						static_cast<const void *>(a), lda, static_cast<const void *>(b), ldb,
						static_cast<const void *>(&beta), static_cast<void *>(c), ldc);
		}
		
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		T *  // M
		pointer(const matrix<T,F,A> &m)
		{
			//return &m(0,0);
			return const_cast<T*>(&m.data().begin()[0]);
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		T *  // trans(M)
		pointer(const matrix_unary2<matrix<T,F,A>,scalar_identity<T> > &m)
		{
			return pointer(m.expression().expression());
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		T *  // trans(M)
		pointer(const matrix_unary2<const matrix<T,F,A>,scalar_identity<T> > &m)
		{
			return pointer(m.expression().expression());
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		T *  // trans(conj(M))
		pointer(const matrix_unary2<matrix_unary1<matrix<T,F,A>,scalar_conj<T> > const ,scalar_identity<T> > &m)
		{
			return pointer(m.expression().expression().expression());
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		T *  // conj(trans(M))
		pointer(const matrix_unary1<matrix_unary2<matrix<T,F,A>,scalar_identity<T> >,scalar_conj<T> > &m)
		{
			return pointer(m.expression().expression().expression());
		}
		
		template<class E>
		MKL_BOOST_UBLAS_INLINE
		MKL_INT
		leading_dimension(const matrix_expression<E> &e, row_major_tag)
		{
			return e().size1();
		}
		template<class E>
		MKL_BOOST_UBLAS_INLINE
		MKL_INT
		leading_dimension(const matrix_expression<E> &e, column_major_tag)
		{
			return e().size2();
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		MKL_INT
		leading_dimension(const matrix<T,F,A> &m, row_major_tag)
		{
			return m.size2();
		}
		template<class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		MKL_INT
		leading_dimension(const matrix<T,F,A> &m, column_major_tag)
		{
			return m.size1();
		}
		
		template<class T>
		MKL_BOOST_UBLAS_INLINE
		CBLAS_ORDER
		storage_layout(T);
		template<>
		MKL_BOOST_UBLAS_INLINE
		CBLAS_ORDER
		storage_layout(row_major_tag) {
			return CblasRowMajor;
		}
		template<>
		MKL_BOOST_UBLAS_INLINE
		CBLAS_ORDER
		storage_layout(column_major_tag) {
			return CblasColMajor;
		}
		
		template<class T>
		MKL_BOOST_UBLAS_INLINE
		int supported_type(T)
		{ return 0; }
		template<>
		MKL_BOOST_UBLAS_INLINE
		int supported_type(float)
		{ return 1; }
		template<>
		MKL_BOOST_UBLAS_INLINE
		int supported_type(double)
		{ return 2; }
		template<>
		MKL_BOOST_UBLAS_INLINE
		int supported_type(std::complex<float>)
		{ return 3; }
		template<>
		MKL_BOOST_UBLAS_INLINE
		int supported_type(std::complex<double>)
		{ return 4; }
		
		template<class E1, class E2, class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		void
		gemm(const  CBLAS_TRANSPOSE transa, const  CBLAS_TRANSPOSE transb,
			 const matrix_expression<E1> &a, const matrix_expression<E2> &b,
			 matrix<T,F,A> &c)
		{
			typedef T value_type;
			typedef typename F::orientation_category orientation_category;
			if(supported_type(value_type()) > 0) {
				MKL_INT m = a().size1();
				MKL_INT k = a().size2();
				MKL_INT n = b().size2();
				MKL_INT lda = leading_dimension(a(),orientation_category());
				MKL_INT ldb = leading_dimension(b(),orientation_category());
				MKL_INT ldc = leading_dimension(c,orientation_category());
				value_type alpha = OneElement(value_type());  // boost/numeric/ublas/detail/concepts.hpp
				value_type beta = ZeroElement(value_type());  // boost/numeric/ublas/detail/concepts.hpp
				CBLAS_ORDER layout = storage_layout(orientation_category());
				T *pa = pointer(a());
				T *pb = pointer(b());
				T *pc = pointer(c);
				//std::cerr << "effettuo il prodotto gemm";
				gemm(layout, transa, transb, m, n, k, alpha, pa, lda, pb, ldb, beta, pc, ldc);
				//std::cerr << "	ok" << std::endl;
			} else {
				c = gemm(a,b);
			}
		}
		template<class E2, class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		void
		//		dcsrmm(const compressed_matrix<T> & a, const matrix_expression<E2> &b,
		dcsrmm(const MATRICI_SPARSE &a, const matrix_expression<E2> &b,
			   matrix<T,F,A> &c)
		{
			char transa = 'n';
			char matdescra[6];
			matdescra[0] = 'g';
			matdescra[1] = 'l';
			matdescra[2] = 'n';
			matdescra[3] = 'c';
			
			typedef T value_type;
			typedef typename F::orientation_category orientation_category;
			MKL_INT m = a.size1();
			MKL_INT k = a.size2();
			MKL_INT n = b().size2();
			MKL_INT ldb = leading_dimension(b(),orientation_category());
			MKL_INT ldc = leading_dimension(c,orientation_category());
			value_type alpha = OneElement(value_type());
			value_type beta = ZeroElement(value_type());
			value_type *value = const_cast<value_type*>(a.value_data().begin());
			value_type *pb = pointer(b());
			value_type *pc = pointer(c);
			MKL_INT rowIndex[m+1];
			int i;
			MKL_INT ri;
			for ( i = 0; i < a.filled1(); ++i) {
				ri = static_cast<MKL_INT> ((a.index1_data()[i]));
				rowIndex[i] = ri;
			}
			for ( ; i < m+1; ++i) {
				ri = static_cast<MKL_INT> (a.nnz());
				rowIndex[i] = ri;
			}
			
			mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, value,(MKL_INT*) &(a.index2_data()[0]),  &(rowIndex[0]), &(rowIndex[1]), pb, &ldb, &beta, pc, &ldc);	
		}
		
		template<class E2, class T, class F, class A>
		MKL_BOOST_UBLAS_INLINE
		void
		dcsrmm_plus(const MATRICI_SPARSE &a, const matrix_expression<E2> &b,
					matrix<T,F,A> &c)
		{
			char transa = 'n';
			char matdescra[6];
			matdescra[0] = 'g';
			matdescra[1] = 'l';
			matdescra[2] = 'n';
			matdescra[3] = 'c';
			
			typedef T value_type;
			typedef typename F::orientation_category orientation_category;
			MKL_INT m = a.size1();
			MKL_INT k = a.size2();
			MKL_INT n = b().size2();
			MKL_INT ldb = leading_dimension(b(),orientation_category());
			MKL_INT ldc = leading_dimension(c,orientation_category());
			value_type alpha = OneElement(value_type());
			value_type beta = OneElement(value_type());
			//value_type beta = ZeroElement(value_type());
			value_type *value = const_cast<value_type*>(a.value_data().begin());
			value_type *pb = pointer(b());
			value_type *pc = pointer(c);
			MKL_INT rowIndex[m+1];
			int i;
			MKL_INT ri;
			for ( i = 0; i < a.filled1(); ++i) {
				ri = static_cast<MKL_INT> ((a.index1_data()[i]));
				rowIndex[i] = ri;
			}
			for ( ; i < m+1; ++i) {
				ri = static_cast<MKL_INT> (a.nnz());
				rowIndex[i] = ri;
			}
			
			mkl_dcsrmm(&transa, &m, &n, &k, &alpha, matdescra, value,(MKL_INT*) &(a.index2_data()[0]),  &(rowIndex[0]), &(rowIndex[1]), pb, &ldb, &beta, pc, &ldc);	
		}
		
	}  // namespace mkl
	
	
    template<class T, class F, class A>
    MKL_BOOST_UBLAS_INLINE
    matrix<T,F,A>    // prod( m1, m2 )
    prod_mkl(const matrix<T,F,A> &m1,
             const matrix<T,F,A> &m2)
    {
        matrix<T,F,A> temporary(m1.size1(),m2.size2());
        mkl::gemm(CblasNoTrans, CblasNoTrans, m1, m2, temporary);
        return temporary;
    }
	template<class T, class F, class A>
	MKL_BOOST_UBLAS_INLINE
    matrix<T,F,A>    // prod( trans(m1), m2 )
    prod_mkl_trans1(const matrix_unary2<const matrix<T,F,A>,scalar_identity<T> > &m1,
         const matrix<T,F,A> &m2)
    {
        matrix<T,F,A> temporary(m1.size1(),m2.size2());
        mkl::gemm(CblasTrans, CblasNoTrans, m1, m2, temporary);
        return temporary;
    }	
	template<class T, class F, class A>
    MKL_BOOST_UBLAS_INLINE
    matrix<T,F,A>    // prod( m1, m2 )
    prod_mkl_sparse(const MATRICI_SPARSE &m1,
					const matrix<T,F,A> &m2)
    {
		
        matrix<T,F,A> temporary(m1.size1(),m2.size2());
		mkl::dcsrmm(m1,m2,temporary);
		return temporary;
    }
	template<class T, class F, class A>
    MKL_BOOST_UBLAS_INLINE
    matrix<T,F,A>    // prod( m2, m1 )
    prod_mkl_sparse_trans(const MATRICI_SPARSE &m1,
						  const matrix<T,F,A> &m2)
    {
 	
		matrix<T,F,A> temporary(m1.size1(),m2.size2());
		mkl::dcsrmm(m1,m2,temporary);
		return trans(temporary);
    }
	template<class T, class F, class A>
    MKL_BOOST_UBLAS_INLINE
    void    // prod( m1, m2 )
    prod_mkl_sparse(const MATRICI_SPARSE &m1,
					const matrix<T,F,A> &m2,
					matrix<T,F,A> &v)
    {
		mkl::dcsrmm_plus(m1,m2,v);
    }
	
/*	template<class T, class F, class A>
    MKL_BOOST_UBLAS_INLINE
    void    // prod( m2, m1 )
    prod_mkl_sparse_trans(const MATRICI_SPARSE &m1,
						  const matrix<T,F,A> &m2,
						  matrix<T,F,A> &v)
    {
		
		matrix<T,F,A> temporary(m1.size1(),m2.size2());
		mkl::dcsrmm(m1,m2,temporary);
		return trans(temporary);
    }
*/	
}}}
#endif  // BOOST_VERSION
#endif  // NDEBUG
#endif  // _MKL_BOOST_UBLAS_MATRIX_PROD_
