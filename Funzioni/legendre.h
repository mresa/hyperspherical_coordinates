/*
 *  legendre.h
 *
 *  Legendre polynomials and derivatives P_l(x)
 *
 *
 */

/*
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

#ifndef LEG_POLI_H
#define LEG_POLI_H

#include <boost/math/special_functions/legendre.hpp>
#include <assert.h>

//P_l(x)
inline long double legendre(const unsigned short & l,
							const long double & x) {
	
	assert (std::fabs(x) <= 1.);                                    //controllo su: ordine l positivo , x tra -1 e 1
	if (l==0) return 1.;													
	else if (l==1) return x;													
	else return boost::math::legendre_p((int)l , x);								
}

//P'_l(x)
inline long double dlegendre(const unsigned short & l,
							 const long double & x) {               // calcola la derivata prima P'_l(x)
	assert (std::fabs(x) < 1.);										// controllo dei limiti: P'_l(x) può divergere agli estremi
	if (l==0) return 0.;											// P'_0 == 0
	else if (l==1) return 1.;										// P'_1(x) == 1
	else {                                                          // calcolo iterativo di P'_l(x)
		return (l*legendre(l-1, x) - l*x*legendre(l,x))/(1-x*x);
	}
}

#endif
