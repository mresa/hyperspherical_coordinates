/*
 *  interfaces_ublas.h.h
 *  Tesi
 *
 *  Created by Marco Resa on 02/08/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef IETL_UBLAS_H
#define IETL_UBLAS_H

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include "tipi_definiti.h"
#include "algebra_matrici.h"
#include <ietl/traits.h>
#include <omp.h> //calcolo parallelo OMP


#include <boost/timer.hpp>

namespace ietl {
	
	template < class T, class Gen> 
    inline void generate_begin(boost::numeric::ublas::vector<T> & c, Gen& gen) {		
		std::generate(c.begin(),c.end(),gen);
	}  
  	
	template < class T, class Gen> 
    inline void regenerate_begin(boost::numeric::ublas::vector<T> & c, Gen& gen, const int & D) {
		int D0 = c.size();
		c.resize(D);
#pragma omp parallel for schedule( guided )
		for (int i = D0; i < D ; ++i) {
			c(i) = abs_tol_0;
		}
	}  

	template < class T, class Gen> 
    inline void generate(boost::numeric::ublas::vector<T> & c, Gen& gen) {		
		c.assign(gen);
	}  
	
	template < class T> 
    inline void clear(boost::numeric::ublas::vector<T>& c) {
		c.clear();
	}  
	
	
/*	template < class T, class S>
	inline T dot(const boost::numeric::ublas::vector<T,S>& x , const boost::numeric::ublas::vector<T,S>& y) {
		return boost::numeric::ublas::inner_prod (boost::numeric::ublas::conj(x), y);
	}
*/	
	template < class T, class S>
	inline T dot(const boost::numeric::ublas::vector<T,S>& x , const boost::numeric::ublas::vector<T,S>& y) {
		return boost::numeric::ublas::inner_prod (x, y);
	}

	template < class T>
	inline typename number_traits<T>::magnitude_type two_norm(const boost::numeric::ublas::vector<T>& x) {
		return boost::numeric::ublas::norm_2(x);
	}
	
	template < class T>
	void copy(const boost::numeric::ublas::vector<T>& x,boost::numeric::ublas::vector<T>& y) {
		y.assign(x);
	}
	
/*	template < class T>
	void mult(const boost::numeric::ublas::matrix<T>& m, 
			  const boost::numeric::ublas::vector<T>& x, 
			  boost::numeric::ublas::vector<T>& y) {
		y=boost::numeric::ublas::prod(m,x);
	}
	
	void mult(const MATRICI_SPARSE & m, 
			  const VETTORI & x, 
			  VETTORI & y) {
		boost::timer t;
		axpy_prod (m,x,y);
		//y=boost::numeric::ublas::prod(m,x);
		std::cerr << "MULT  T = " << t.elapsed() << std::endl;
	}
*/	
	void mult(const V_V_MATRICI_SPARSE & M,
			  const VETTORI & x,
			  VETTORI & y) {
		assert(M.size() > 0);
		y.assign(v_mult(M,x));
	}
	
}

#endif
