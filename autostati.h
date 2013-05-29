/*
 *  autostati.h
 *  Tesi
 *
 *  Created by Marco Resa on 28/11/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef AUTOSTATI_H
#define AUTOSTATI_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/vector.hpp>
#define BOOST_NUMERIC_BINDINGS_USE_CLAPACK
#include <boost/numeric/bindings/lapack/syevd.hpp>
//#include <boost/numeric/bindings/lapack/syev.hpp>
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#undef  BOOST_NUMERIC_BINDINGS_USE_CLAPACK


void autostati(const boost::numeric::ublas::symmetric_matrix<tipo_matrici>& A, boost::numeric::ublas::matrix<tipo_matrici>& Q, boost::numeric::ublas::vector<tipo_matrici>& d)
{
    namespace ublas  = boost::numeric::ublas;
    namespace lapack = boost::numeric::bindings::lapack;
	
    ublas::vector<tipo_matrici> cd(A.size1());
    ublas::matrix<tipo_matrici, ublas::column_major> CQ(A.size1(), A.size2());
    int info;
	
    for (std::size_t i = 0; i < A.size1(); ++i) {
        for (std::size_t j = 0; j <= i; ++j) {
            CQ(i, j) = A(i, j);
        }
    }
	
    info = lapack::syevd('V', 'L', CQ, cd, lapack::optimal_workspace());
    BOOST_UBLAS_CHECK(info == 0, ublas::internal_logic());
	
    Q = CQ;
    d.swap(cd);
}


#endif
