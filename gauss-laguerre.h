/*
 *  gauss.h
 *
 *	Tesi
 *  Abramowitz : 25.4.45 integrazione con peso e^-x tra 0 e inf
 *  
 *	La classe Gauss_Laguerre dato l'ordine n, alpha costruisce 2 vettori composti dagli zeri di L_n^alpha(x) e dai pesi della formula di quadratura di gauss-laguerre
 *
 *  Viene creata con 2 metodi:
 *  G(n)  (alpha == 0)
 *  G(n,alpha)
 *
 *
 *	Fornisce le funzioni pubbliche
 *	N() che restituisce l'ordine n
 *	ALPHA() che restituisce l'ordine alpha
 *	X(n) restituisce lo zero n-esimo
 *	W(n) restituisce il peso n-esimo
 *
 *	Per l'algoritmo di newton richiede il parametro long double err, valore di stop del decremento
 *
 *
 */

#ifndef GAUSS_LAG_H
#define GAUSS_LAG_H
#include "err.h"
#include "laguerre.h"
#include <assert.h>
// 25.4.45

class Gauss_Laguerre {
	unsigned short n;						//ordine laguerre
	long double* x;							//griglia degli zeri
	long double* w;							//griglia dei pesi
	unsigned short alpha;					//laguerre generalizzato
	
	static inline bool errore_gauss(){
		std::cerr << "errore in gauss-laguerre.h" << std::endl;
		std::exit(1);
		return 1;
	}

	static inline bool errore_zeri(){
		std::cerr << "errore negli zeri - gauss - laguerre" << std::endl;
		std::exit(1);
		return 1;
	}
	
	static inline void zero(const unsigned short & n, 
							const unsigned short & alpha, 
							long double* zeri) {
		if ( n<1 ) {													// controllo dei limiti: L_n(x) ha n zeri; in particolare L_0 non ne ha
			errore_gauss();
		}
		else {
			
			/*
			 *	 Guess iniziale per gli zeri: referenze "Numerical recipes in fotran77 ...", Cambridge U Press, 1999 oppure gen_laguerre_rule.C su internet
			 */
			
			for (short i = 0 ; i < n ; ++i) {  //NB non posso parallelizzarlo, ognuno dipende dal precedente
				if ( i == 0 ) {
					zeri[i] = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) / 
					( 1.0 + 2.4 * n + 1.8 * alpha );
					
				}
				else if ( i == 1 ) {
					zeri[i] = zeri[i-1] + ( 15.0 + 6.25 * alpha ) / 
					( 1.0 + 0.9 * alpha + 2.5 * n );
				}
				else {
					long double r1 = ( 1.0 + 2.55 * ( long double ) ( i - 1 ) ) 
					/ ( 1.9 * (long double ) ( i - 1 ) );
					
					long double r2 = 1.26 * ( long double ) ( i - 1 ) * alpha / 
					( 1.0 + 3.5 * ( long double ) ( i - 1 ) );
					
					long double ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );
					
					zeri[i] = zeri[i-1] + ratio * ( zeri[i-1] - zeri[i-2] );
				}
				
				/*
				 *	Algoritmo newton per il calcolo degli zeri: si ferma per incrementi minori di err
				 */
				
				while (std::fabs(laguerre(n, alpha, zeri[i])/dlaguerre(n, alpha, zeri[i])) > (err*zeri[i]) && std::fabs(laguerre(n, alpha, zeri[i])/dlaguerre(n, alpha, zeri[i])) > err){
						zeri[i] -= (laguerre(n, alpha, zeri[i])/dlaguerre(n, alpha, zeri[i]));
				}
				
				
				/* 
				 *  controllo dei limiti per lo zero trovato
				 */

				if (controllo_zeri_laguerre(n, alpha, zeri[i], i+1)) ;
				else errore_zeri();
			}
		}
	}
	
	static inline void zero(const unsigned short & n, 
							long double* zeri) {
		zero(n, 0, zeri);
	}
	
	
public:
	
	
	Gauss_Laguerre (const unsigned short & N, 
					const unsigned short & ALPHA ) {
		if (N <= 0) errore_gauss();
		else {
			n = N;
			alpha = ALPHA;
			x = new long double[n];
			w = new long double[n];
			zero(n , alpha, x);

//			Calcolo dei pesi W_i
			
			for ( short i=0 ; i < n ; ++i) {

				/*
				 *		la funzione boost::math::rising_factorial( n + 1 , alpha ) == (n + alpha)! / n! == (n + 1)(n + 2) … (n + alpha)
				 *		la funzione boost::math::tgamma_delta_ratio = 1/ boost::math::rising_factorial
				 *
				 *
				 */
				
//				w[i] = boost::math::rising_factorial( n + 1 , alpha ) * (x[i]  /  pow( ( (n+1)*laguerre(n+1, alpha, x[i]) ) ,2));
				w[i] = x[i]  / ( pow( ( (n+1)*laguerre(n+1, alpha, x[i]) ) ,2) * boost::math::tgamma_delta_ratio( n + 1 , alpha ));
			}
		}	
	}	
	
	
	Gauss_Laguerre ( const unsigned short & N ) {
		if (N <= 0) errore_gauss();
		else {
			n = N;
			alpha = 0;
			x = new long double[n];
			w = new long double[n];
			zero(n , x);

			//			Calcolo dei pesi W_i
			
			for (short i=0 ; i<n ; ++i) {
				w[i] = x[i] / ( pow( (n+1)*laguerre(n+1,x[i]) , 2 ) );
			}
		}	
	}		

	
	//distruttore

	~Gauss_Laguerre () {
		delete[] x;
		delete[] w;
	}
	
	inline long double X(const unsigned short & i) const {
		assert( i < n );
		return x[i];
	}
	inline long double W(const unsigned short & i) const {
		assert( i < n );
		return w[i];
	}
	inline unsigned short N() const {return n;}
	inline long double ALPHA() const {return alpha;}
};


#endif
