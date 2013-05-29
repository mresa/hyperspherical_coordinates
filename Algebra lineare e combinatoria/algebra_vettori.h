/*
 *  algebra_vettori.h
 *  Tesi
 *
 *  Created by Marco Resa on 05/08/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef ALGEBRA_VECTOR_H
#define ALGEBRA_VECTOR_H

#include <iostream>
#include <vector>
#include <omp.h> //calcolo parallelo OMP

template <class T>
inline bool delta_vec(const std::vector <T> & l1, 
					  const std::vector <T> & l2){
	assert (l1.size() == l2.size());
	bool temp = true;
	if (l1.size() != l2.size()) {
		temp = false;
	}
	else{
		for (short i = 0; temp *(i < l1.size()) ; ++i){
			temp *= l1[i]==l2[i];
		}		
	}
	return temp;
}

template <class T>
inline bool delta_vec_tranne(const std::vector <T> & l1, 
							 const std::vector <T> & l2,
							 const unsigned short & p1,
							 const unsigned short & p2){
	assert (l1.size() == l2.size());
	bool temp = true;
	if (l1.size() != l2.size()) {
		temp = false;
	}
	else {
		for (short i = 0; temp *(i < l1.size()) ; ++i){
			if ((i != p1)&&(i != p2)) {
				temp *= l1[i]==l2[i];
			}
		}		
	}

	return temp;
}

template <class T>
inline bool delta_vec_fino(const std::vector <T> & l1, 
						   const std::vector <T> & l2,
						   const unsigned short & p1){
	assert (l1.size() == l2.size());
	bool temp = true;
	if (l1.size() != l2.size() || p1 > l1.size()) {
		temp = false;
	}
	else {
		for (short i = 0; temp *(i < p1) ; ++i){
			temp *= l1[i]==l2[i];
		}		
	}
	
	return temp;
}

/*template <class T>
inline T quad_sum (std::vector <T> coefficienti) {
	T s = 0;
	for(std::vector <T> ::iterator i = coefficienti.begin(); i != coefficienti.end() ; ++i){		
		s += std::pow((*i),2);
	}
	return s;
}
*/
template <class T>
inline T quad_sum (std::vector <T> coefficienti) {
	T s = 0;
#pragma omp parallel for reduction(+:s) schedule( guided )
	for(short i = 0; i < coefficienti.size() ; ++i){		
		s += (coefficienti[i])*(coefficienti[i]);
	}
	return s;
}
 
template <class T>
inline T linear_sum (std::vector <T> coefficienti) {
	T s = 0;
#pragma omp parallel for reduction(+:s) schedule( guided )
	for(short i = 0; i < coefficienti.size() ; ++i){		
		s += (coefficienti[i]);
	}
	return s;
}

inline unsigned short linear_sum_bool (std::vector <bool> coefficienti) {
	unsigned short s = 0;
#pragma omp parallel for reduction(+:s) schedule( guided )
	for(short i = 0; i < coefficienti.size() ; ++i){		
		s += (coefficienti[i]);
	}
	return s;
}

#endif /* ALGEBRA_VECTOR_H */

