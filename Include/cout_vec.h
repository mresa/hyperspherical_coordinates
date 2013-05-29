/*	
 *  cout_vec.h
 *  Tesi
 *	
 *	Funzione di supporto per il cout dei vettori tramite cout_vec(vettore)
 *	Stampa a schermo un vettore
 *	
 *
 *  Created by Marco Resa on 25/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef COUT_VEC_H
#define COUT_VEC_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

template <class BidirectionalIterator>
inline void cout_vec_support(const BidirectionalIterator & primo, const BidirectionalIterator & ultimo) {		//boleano, vuole l'iteratore corrispondente al primo e all'ultimo elemento
	BidirectionalIterator corrente = primo;
	if(primo != ultimo){
		do{
			std::cout << *corrente << " ";
		}while (++corrente != ultimo);
	}
	else;
}

template <class BidirectionalIterator>
inline void cerr_vec_support(const BidirectionalIterator & primo, const BidirectionalIterator & ultimo) {		//boleano, vuole l'iteratore corrispondente al primo e all'ultimo elemento
	BidirectionalIterator corrente = primo;
	if(primo != ultimo){
		do{
			std::cerr << *corrente << " ";
		}while (++corrente != ultimo);
	}
	else;
}

template <class T>
inline void cout_vec(const std::vector<T> & V) {
	if(V.size()!=0){
		std::cout << " [ ";
		cout_vec_support(V.begin(), V.end());
		std::cout << "] ";			
	}
	else;
}

template <class T>
inline void cerr_vec(const std::vector<T> & V) {
	if(V.size()!=0){
		std::cerr << " [ ";
		cerr_vec_support(V.begin(), V.end());
		std::cerr << "] ";			
	}
	else;
}

#endif
