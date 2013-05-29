/*
 *  err.h
 *  Tesi
 *	Errore utilizzato nella ricerca degli zeri per le quadrature; migliora la precisione a prezzo di tempo macchina
 *  Created by Marco Resa on 16/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef ERR_H
#define ERR_H 
#include <limits>


const long double err = 100*std::numeric_limits<long double>::epsilon(); 
// 1.1e-19 x 100 sulla mia macchina, 2 ordini di grandezza l'epsilon macchina
const long double err_coordinate = 1000*std::numeric_limits<long double>::epsilon();
// 1.1e-19 x 1000 sulla mia macchina, 2 ordini di grandezza l'epsilon macchina
																					
const long double err_matrici = 100000*err_coordinate;
//const long double err_matrici = 1000000000*err_coordinate;
const double abs_tol_0 = 10*std::pow(std::numeric_limits<double>::epsilon(),2./3); 
//3.6e-10
const long double errore_matrici = 1.1*abs_tol_0;
const long double errore_moltiplicazione_matrici = 10000000*std::numeric_limits<long double>::epsilon(); 
// 1.1e-12

/*
const long double err = 100*std::numeric_limits<long double>::epsilon(); 
const long double err_coordinate = 1000*std::numeric_limits<long double>::epsilon();

const long double err_matrici = err_coordinate;
const double abs_tol_0 = 10*std::pow(std::numeric_limits<double>::epsilon(),2./3); 
const long double errore_matrici = err_coordinate;
const long double errore_moltiplicazione_matrici = err_coordinate; 
// 1.1e-12
*/

#endif

