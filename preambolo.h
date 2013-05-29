/*
 *	preambolo.h
 *
 *	Copyright (C) 2013 by Marco Resa <https://github.com/mresa>
 *
 *	This program is free software: you can redistribute it and/or modify
 *	it under the terms of the GNU General Public License as published by
 *	the Free Software Foundation, either version 3 of the License, or
 *	(at your option) any later version.
 *
 *	This program is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *	GNU General Public License for more details.
 *
 *	You should have received a copy of the GNU General Public License
 *	along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


//Parametri di ottimizzazione XCode


#ifndef PREAMBOLO_H
#define PREAMBOLO_H

//#if defined(XCODE_DEBUG_CONFIGURATION_DEBUG)
//#define DEBUG 1
//#else
//#define DEBUG 1
//#define NDEBUG 1
//#define BOOST_UBLAS_NDEBUG 1
//#define OMP_NUM_THREADS 8
//#define MKL_NUM_THREADS 8
//#endif

#include <cassert>
// input output
#include <iostream>
// scrittura file
#include <fstream>
//mathematica
#include <cmath>
//gestione limiti computazionali
#include <limits>
//vettori
#include <vector>
//stampa un vettore su cout
#include "cout_vec.h"
//libreria stringhe
#include "stringhe.h"
//libreria di algebra e algebra lineare

// commentato a posteriori
//#include <mkl_service.h>

//#include <mkl.h>

#include "algebra.h"

#include "algebra_matrici.h"

#include "algebra_vettori.h"

//funzione di bessel
#include <gsl/gsl_sf_bessel.h>
//funzione di laguerre
#include <boost/math/special_functions/laguerre.hpp>
//gamma di eulero tgamma(n+1)=n!
#include <boost/math/special_functions/gamma.hpp>
//fattoriali
#include <boost/math/special_functions/factorials.hpp>

#include <sys/stat.h>

#endif
