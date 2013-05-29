/*
 *  costanti_fisiche.h
 *  Tesi
 *
 *  Created by Marco Resa on 01/07/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */


#ifndef COSTANTI_FISICHE_H
#define COSTANTI_FISICHE_H

#include <boost/math/constants/constants.hpp>
//pi greco in long double
long double pi = boost::math::constants::pi <long double> ();
long double planck_h = 4.13566733e-15 ;
long double h_bar = 6.58211899e-16 ;
long double struttura_fine = 1/137.035999679 ;
long double h_bar_c = 197.3269631;
long double e_quadro = 1.44;

//massa neutronica
//long double m_n = 939.565 = h_bar_c*h_bar_c/41.44
long double m_n = h_bar_c*h_bar_c/41.44;

//massa protonica
//long double m_p = 938.272 = h_bar_c*h_bar_c/41.50
long double m_p = h_bar_c*h_bar_c/41.50;

//938.94213567 = massa di riferimento
long double m_ref = h_bar_c*h_bar_c/41.47;
//m_ref = 2*m_n*m_p/(m_n + m_p) = h_bar_c*h_bar_c/41.47

// conversione MeV -> u unità atomiche
long double conv_u =  1.00727646677 / 938.272013;
long double beta_riferimento = 3.;
long double m_epsilon = (m_n - m_p)/(m_n + m_p); //0.0006923

long double m_A = m_ref*(m_n + m_p)/(2*m_n*m_p) - 1;
long double m_B = m_ref*(m_n - m_p)/(2*m_n*m_p);

#endif
