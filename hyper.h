/*
 *  hyper.h
 *  Tesi
 *	Crea l'armonica ipersferica dati i numeri quantici.
 *
 *	Creazione tramite vettori
 *
 *			v_n, v_l, v_m
 *
 *			v_n
 * Creazione come potential basis
 *			N , n, l, m
 *			N , n
 * 
 *
 *
 *
 *
 *  Created by Marco Resa on 05/11/08.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef HYPER_H
#define HYPER_H
#include <assert.h>
#include <vector>
#include "tipi_definiti.h"
/*	
 *	porta in dotazione std::complex<calculated-result-type> spherical_harmonic(unsigned n, int m, T1 theta, T2 phi);
 *	θ in [0, π], φ in [0,2π)
 *
 */
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include "algebra.h"
#include "jacobi.h"

//armonicha ipersferica di ordine 0 per N dato
inline long double hyper0 (const unsigned short & N){
	return sqrt((long double) boost::math::tgamma((long double)3.*N/2.)/(2*pow(pi,(long double)3.*N/2.)));
}

inline long double potential_basis_normalizzazione (const unsigned short& N,
													const unsigned short& n){
	long double nu = 2. * n + 3.* N / 2. - 1.;
	long double normalizzazione = sqrt(
									   2.* nu * 
									   boost::math::tgamma_delta_ratio(n + 1, 0.5)
									   /
									   boost::math::tgamma_delta_ratio(nu - n - 0.5, 0.5));
	normalizzazione *= 2. * sqrt(pi) / (sqrt(pow((long double)2., 3 * N )) * hyper0(N-1) );
	
	return normalizzazione;
}

inline long double potential_basis_normalizzazione (const unsigned short& N,
													const unsigned short& n,
													const unsigned short& l){
	long double nu = 2. * n + l + 3.* N / 2. - 1.;
	long double normalizzazione = sqrt(
									   2.* nu * 
									   boost::math::tgamma_delta_ratio(n + 1., l + 1./2.)/
									   boost::math::tgamma_delta_ratio(nu - n - l - 1./2., l + 1./2.));
	return normalizzazione;
}

inline long double tre_corpi_normalizzazione (const unsigned short& n,
											  const unsigned short& l1,
											  const unsigned short& l2){
	unsigned short K = 2. * n + l1 + l2 ;
	long double normalizzazione = sqrt(
									   2.* (K+2) * 
									   boost::math::tgamma_ratio(K + 2 - n , n + l1+ 3./2.) * 
									   boost::math::tgamma_ratio(n + 1, n + l2 + 3./2.));
	return normalizzazione;
}

//classe per la gestione delle armoniche ipersferiche
class hyper {

private:
	
	unsigned short nn;									// N = particelle -1
	std::vector< unsigned short > n_quantici;			// sono N : 0 , [n_2 … n_N] : n_1 == 0
	std::vector< unsigned short > l_quantici;			// sono N : [l_1 … l_N]			
	std::vector< short > m_quantici;					// sono N : [m_1 … m_N]
	std::vector< unsigned short > k_quantici;			// sono N
	std::vector< long double > nu_quantici;				// sono N
	std::vector< long double > normalizzazione;			//normalizzazione per l'indice i
	
	
		
	//controllo : n_1 == 0
	static void controlla_n(const std::vector< unsigned short >&  NN ) {
		assert (NN[0] == 0);
	}

	// inizializza il gruppo n
	static void init_n (const std::vector< unsigned short >&  NN , 
						std::vector< unsigned short > *  nn ) {
		(*nn).resize(NN.size() + 1);
		(*nn)[0]=0.;
		for (int i = 1 ; i < (*nn).size() ; ++i) {
			(*nn)[i] = NN[i-1];
		}
	}
	
	// inizializza il gruppo K = sum(2n  + l)
	static void init_K (const std::vector< unsigned short >&  NN , 
						const std::vector< unsigned short >&  LL , 
						std::vector< unsigned short > * kk) {
		(*kk).resize(NN.size());
		(*kk)[0]=LL[0];
		for (int i = 1 ; i < (*kk).size() ; ++i) {
			(*kk)[i] = (*kk)[i-1] + 2*NN[i] + LL[i];
		}
	}
	
	
	// inizializza il gruppo nu_(i+1) = k_(i+1) + 3i/2 - 1 
	static void init_nu (const std::vector< unsigned short > & KK , 
						 std::vector< long double > * nu ) {
		(*nu).resize(KK.size());
		for (int i = 0 ; i < (*nu).size() ; ++i) {
			(*nu)[i] = KK[i] + 3. * (i + 1) /2. - 1.;
		}
	}
	
	// inizializza la normalizzazione
	static void init_normalizzazione (const std::vector< unsigned short > & n,
									  const std::vector< unsigned short > & l,
									  const std::vector< unsigned short > & k,
									  const std::vector< long double > & nu,
									  std::vector< long double > * normalizzazione) {
		(*normalizzazione).push_back(1.);
		
		for (int i = 1 ; i < n.size() ; ++i) {
			//giusto un'altra versione
/*			(*normalizzazione).push_back(sqrt(
											  2.*nu[i]* pochhammer(nu[i] - n[i] - l[i] - 1./2., l[i] + 1./2.)/
											  pochhammer(n[i] + 1., l[i] + 1./2.)));
*/		
			(*normalizzazione).push_back(sqrt(
											  2.*nu[i]*
											  boost::math::tgamma_delta_ratio(n[i] + 1., l[i] + 1./2.)/
											  boost::math::tgamma_delta_ratio(nu[i] - n[i] - l[i] - 1./2., l[i] + 1./2.)));
		}
	}

	// Funzione per il calcolo dell'armonica sferica all'indice i, dati gli angoli
	static inline complesso armonica_sferica (const unsigned short & i,
									   const std::vector< unsigned short >& l,
									   const std::vector< short >& m,
									   const long double & teta, 
									   const long double & phi ) {
		return boost::math::spherical_harmonic(l[i], m[i], teta, phi);
	}

	
	// Funzione per il calcolo della funzione di ripelle
	static inline long double ripelle (const unsigned short& i,
									   const std::vector< unsigned short >& n,
									   const std::vector< unsigned short >& l,
									   const std::vector< unsigned short >& k,
									   const std::vector< long double >& nu,
									   const long double& FI ) {
		if(i==0) return 1. ;
		return pow( cos( FI ) , l[i] ) * pow( sin( FI ) , k[i-1] ) * jacobi( n[i] , nu[i-1] , l[i] + 1./2. , cos(2.*FI) ) ; 
	}

	static inline long double ripelle_z (const unsigned short& i,
										 const std::vector< unsigned short >& n,
										 const std::vector< unsigned short >& l,
										 const std::vector< unsigned short >& k,
										 const std::vector< long double >& nu,
										 const long double& z ) {
		if(i==0) return 1. ;
		return pow( sqrt((1.+z)/2.) , l[i] ) * pow( sqrt((1.-z)/2.) , k[i-1] ) * jacobi( n[i] , nu[i-1] , l[i] + 1./2. , z ) ; 
	}
	
public:
	
	
	//crea a partire dai numeri quantici
	hyper(const std::vector< unsigned short >&  N_QUANTICI , 
		  const std::vector< unsigned short >&  L_QUANTICI , 
		  const std::vector< short >&  M_QUANTICI ) {
		assert (N_QUANTICI.size() >= 0);
	
		nn = N_QUANTICI.size() + 1;
		
		assert (L_QUANTICI.size() == nn);
		assert (M_QUANTICI.size() == nn);
		
		init_n( N_QUANTICI, &n_quantici);
		l_quantici = L_QUANTICI;
		m_quantici = M_QUANTICI;
		init_K( n_quantici , l_quantici , &k_quantici );
		init_nu( k_quantici , &nu_quantici );
		init_normalizzazione( n_quantici, l_quantici, k_quantici, nu_quantici, &normalizzazione );
		controlla_n (n_quantici);
	}	
	
	//L_i == 0
	hyper ( const std::vector< unsigned short >&  N_QUANTICI ) {
		assert (N_QUANTICI.size() >= 0);
		
		nn = N_QUANTICI.size() + 1;
		
		init_n( N_QUANTICI, &n_quantici);
		l_quantici.assign( nn, 0 );
		m_quantici.assign( nn, 0 );
		init_K( n_quantici , l_quantici , &k_quantici );
		init_nu( k_quantici , &nu_quantici );
		init_normalizzazione( n_quantici, l_quantici, k_quantici, nu_quantici, &normalizzazione );
		controlla_n (n_quantici);
	}	
	
	//crea la potential basis per l'ennesimo elemento
	hyper(const unsigned short& N, 
		  const unsigned short&  N_QUANTICO , 
		  const unsigned short&  L_QUANTICO , 
		  const short&  M_QUANTICO ) {

		assert (N > 0);
		assert (std::abs(M_QUANTICO) <= L_QUANTICO);
		
		nn = N;
		
		n_quantici.assign( nn-1, 0 );
		n_quantici.push_back (N_QUANTICO);
		
		l_quantici.assign( nn-1, 0 );
		l_quantici.push_back (L_QUANTICO);
		
		m_quantici.assign( nn-1, 0 );
		m_quantici.push_back (M_QUANTICO);
		
		init_K( n_quantici , l_quantici , &k_quantici );
		init_nu( k_quantici , &nu_quantici );
		init_normalizzazione( n_quantici, l_quantici, k_quantici, nu_quantici, &normalizzazione );
		controlla_n (n_quantici);
	}	
	
	//crea n_N != 0 ; altri tutti == 0
	hyper( const unsigned short& N, const unsigned short&  N_QUANTICO) {
		assert (N > 0);
		nn = N;
		
		n_quantici.assign( nn-1, 0 );
		n_quantici.push_back (N_QUANTICO);
		
		
		l_quantici.assign( nn, 0 );
		
		m_quantici.assign( nn, 0 );
		
 		init_K( n_quantici , l_quantici , &k_quantici );
		init_nu( k_quantici , &nu_quantici );
		init_normalizzazione( n_quantici, l_quantici, k_quantici, nu_quantici, &normalizzazione );
		controlla_n (n_quantici);
	}	
	
	//distruttore
	~hyper () {
		std::vector< unsigned short >().swap(n_quantici);
		std::vector< unsigned short >().swap(l_quantici);
		std::vector< short >().swap(m_quantici);
		std::vector< unsigned short >().swap(k_quantici);
		std::vector< long double >().swap(nu_quantici);
	}
	
	
	//accesso a N
	unsigned short N() const {return nn;}

	//accesso a K
	unsigned short K() const {return k_quantici[nn - 1];}

	//accesso a [n]
	std::vector< unsigned short > n() const {return n_quantici;}
	unsigned short n( const unsigned short& i) const {
		assert (i < nn);
		return n_quantici[i];
	}

	//accesso a [l]
	std::vector< unsigned short > l() const {return l_quantici;}
	unsigned short l(const unsigned short& i) const {
		assert (i < nn);
		return l_quantici[i];
	}

	//accesso a [m]
	std::vector< short > m() const {return m_quantici;}
	short m(const unsigned short& i) const {
		assert (i < nn);
		return m_quantici[i];
	}

	//accesso a [k]
	std::vector< unsigned short > k() const {return k_quantici;}
	unsigned short k(const unsigned short& i) const {
		assert (i < nn);
		return k_quantici[i];
	}

	//accesso a [nu]
	std::vector< long double > nu() const {return nu_quantici;}
	long double nu(const unsigned short& i) const {
		assert (i < nn);
		return nu_quantici[i];
	}
	
	// alpha : utile per i processi di integrazione, è alpha per la funzione di jacobi
	long double alpha(const unsigned short& i) const {
		assert (i > 0);
		assert (i < nn);
		return nu_quantici[i-1];
	}

	// beta : utile per i processi di integrazione, è beta per la funzione di jacobi
	long double beta(const unsigned short& i) const {
		assert (i < nn);
		return l_quantici[i] + 1./2.;
	}


	//accesso alla normalizzazione per indice i
	std::vector< long double > NORM() const {return normalizzazione;}
	long double NORM( const unsigned short& i ) const {
		assert (i < nn);
		return normalizzazione[i];
	}

	// armonica sferica per indice i, dati gli angoli
	complesso HS (const unsigned short& i, 
				  const long double& teta, 
				  const long double& phi ) {
		assert (i < nn);
		return armonica_sferica( i, l_quantici, m_quantici, teta, phi );
	}
	
	
	// armonica sferica, dati tutti angoli
	complesso HS(const std::vector< long double >& teta, 
				 const std::vector< long double >& phi ) {
		assert (teta.size() == nn);
		assert (phi.size() == nn);
		complesso temp(1.,0.);
		for (short i = 0; i < l_quantici.size() ; ++i) {
			temp *= armonica_sferica(i, l_quantici, m_quantici, teta[i], phi[i]);
		}
		return temp;
	}
	
	// funzione di ripelle normalizzata, indice i, dato l'angolo
	long double RIPELLE(const unsigned short& i,
						const long double& FI ) const {
		assert (i < nn);
		return (normalizzazione[i]
				*
				ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, FI));
	}

	// funzione di ripelle normalizzata, dati tutti gli angoli
	long double RIPELLE(const std::vector< long double >& FI ) const {
		if (FI.size() != nn - 1){
			std::cerr << "FI.size()  " << FI.size() << "  " << nn - 1 << std::endl;
			cerr_vec(FI);
		}
		assert (FI.size() == nn - 1);
		std::vector< long double > fi(FI);
		fi.insert(fi.begin(), 0.); // il primo angolo è sempre zero
		long double temp = 1.;
		for (short i = 0; i < fi.size() && temp != 0. ; ++i) {
			temp *= normalizzazione[i]*ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, fi[i]);
		}
		return temp;
	}
	
	// funzione di ripelle normalizzata, indice i, dato l'angolo
	long double RIPELLE_z(const unsigned short& i,
						  const long double& z ) const {
		assert (i < nn);
		return (normalizzazione[i]
				*
				ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, z));
	}
	
	// funzione di ripelle normalizzata, dati tutti gli angoli
	long double RIPELLE_z(const std::vector< long double >& Z ) const {
		if (Z.size() != nn - 1){
			std::cerr << "Z.size()  " << Z.size() << "  " << nn - 1 << std::endl;
			cerr_vec(Z);
		}
		assert (Z.size() == nn - 1);
		std::vector< long double > z(Z);
		z.insert(z.begin(), 1.); // il primo angolo è sempre zero
		long double temp = 1.;
		for (short i = 0; i < z.size() && temp != 0. ; ++i) {
			temp *= normalizzazione[i]*ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, z[i]);
		}
		return temp;
	}
	// armonica ipersferica, indice i, dati gli angoli
	complesso HHS(const unsigned short& i, 
				  const long double& teta, 
				  const long double& phi,
				  const long double& FI ) {
		assert (i < nn);
		return (armonica_sferica( i, l_quantici, m_quantici, teta, phi ) 
				*
				normalizzazione[i]
				* 
				ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, FI));
	}

	// armonica ipersferica, dati tutti gli angoli
	complesso HHS(const std::vector< long double >& teta, 
				  const std::vector< long double >& phi,
				  const std::vector< long double >& FI)  {
		assert (teta.size() == nn);
		assert (phi.size() == nn);
		assert (FI.size() == nn - 1);		
		std::vector< long double > fi(FI);
		fi.insert(fi.begin(), 0.); // il primo angolo è sempre zero
		complesso temp(1.,0.);
		for (short i = 0; i < fi.size() ; ++i) {
			temp *= armonica_sferica( i, l_quantici, m_quantici, teta[i], phi[i] ) * normalizzazione[i] * ripelle( i, n_quantici, l_quantici, k_quantici, nu_quantici, fi[i]);
		}
		return temp;
	}
		
};

#endif
