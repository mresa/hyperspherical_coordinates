/*
 *  combinazioni.h
 *  Tesi
 *
 *	combinazioni : 
 *	
 *	genera tutti le possibili combinazioni ordinate di n numeri interi positivi a somma fissa K
 *
 *	dati 2 iteratori iniziale e finale, è vero finché esistono altre combinazioni tali che la somma degli elementi sia fissa
 *	si comincia con un vettore fatto di n-1 zeri e K in fondo
 *	
 *	
 *	genera_m : 
 *
 *	genera tutte le possibili combinazioni di n numeri interi che rispettano la regola |m[i]|<=l[i], ovvero genera gli stati possibili di n corpi con l fisso, degeneri su m
 *	
 *	genera_l : 
 *
 *	genera tutte le possibili combinazioni di n numeri interi che rispettano la regola 0 <= l[i]| <= l_max[i]
 *	
 *	dati 2 iteratori iniziale e finale, dato un iteratore finale di un vettore che contiene gli l, è vero finché esistono altre combinazioni degli m per quei valori di l
 *	
 *  Created by Marco Resa on 17/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef COMB_H
#define COMB_H

#include <iterator>
#include <algorithm>


template <class BidirectionalIterator>

inline bool combinazioni(BidirectionalIterator primo, 
				  BidirectionalIterator ultimo) {		//boleano, vuole l'iteratore corrispondente al primo e all'ultimo elemento
	BidirectionalIterator corrente = ultimo;			//il corrente è l'ultimo
	while (corrente != primo && *(--corrente) == 0) {	//finché il corrente è diverso dal primo, itero verso il primo finché non trovo un valore diverso da 0
	}
	if (corrente == primo) {
		if (primo != ultimo && *primo != 0)		//se sono arrivato al primo, se è diverso dall'ultimo (come posizione, non come valore) e se è diverso da zero
			std::iter_swap(--ultimo, primo);	// scambia il valore del penultimo (ora ultimo == penultimo) con il valore del primo e ritorna falso
		return false;
	}
	--(*corrente);						//diminuisce il valore di corrente
	std::iter_swap(--ultimo, corrente); //scambia il penultimo con il corrente
	++(*(--corrente));					//aumenta il valore del penultimo
	return true;
}


template <class BidirectionalIterator>

inline bool genera_m (BidirectionalIterator primo,
			   BidirectionalIterator ultimo,
			   BidirectionalIterator valori) {
	if (ultimo == primo) {
		return false;
	}
	do {
		if (*(--ultimo) > -(*(--valori))) {
 			--(*ultimo);
			return true;
		}
		(*ultimo) = (*valori);
	} while (ultimo != primo);
	return false;
}

template <class BidirectionalIterator>

inline bool genera_m_dub (BidirectionalIterator primo,
						  BidirectionalIterator ultimo,
						  BidirectionalIterator valori) {
	if (ultimo == primo) {
		return false;
	}
	do {
		if (*(--ultimo) > -(*(--valori))) {
 			(*ultimo)-= 2;
			return true;
		}
		(*ultimo) = (*valori);
	} while (ultimo != primo);
	return false;
}

template <class BidirectionalIterator>

inline bool genera_l (BidirectionalIterator primo,
					  BidirectionalIterator ultimo,
					  BidirectionalIterator valori) {
	if (ultimo == primo) {
		return false;
	}
	do {
		if (*(--ultimo) > 0) {
 			--(*ultimo);
			return true;
		}
		(*ultimo) = (*(--valori));
	} while (ultimo != primo);
	return false;
}


#endif /* COMB_H */
