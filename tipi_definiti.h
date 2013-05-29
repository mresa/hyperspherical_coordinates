/*
 *  tipi_definiti.h
 *  Tesi
 *
 *  Created by Marco Resa on 04/08/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TIPI_DEFINITI_H
#define TIPI_DEFINITI_H
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <vector>

typedef std::complex< long double> complesso;

typedef std::vector <short>				vettore_short;
typedef std::vector <unsigned short>	vettore_ushort;
typedef std::vector <long double>		vettore_ldouble;
typedef std::vector <double>			vettore_double;
typedef std::vector <bool>				vettore_bool;

typedef std::vector <vettore_short>		v_vettore_short;
typedef std::vector <vettore_ushort>	v_vettore_ushort;
typedef std::vector <vettore_ldouble>	v_vettore_ldouble;

typedef			double		tipo_matrici;
//typedef			float		tipo_matrici;
//si dichiarano sempre tutti i restanti

//simmetriche
typedef boost::numeric::ublas::symmetric_matrix		<tipo_matrici, boost::numeric::ublas::upper >		MATRICI_SIMMETRICHE_UP;
typedef std::vector									<MATRICI_SIMMETRICHE_UP>							V_MATRICI_SIMMETRICHE_UP;
typedef std::vector									<V_MATRICI_SIMMETRICHE_UP>							V_V_MATRICI_SIMMETRICHE_UP;

//simmetriche nulle in diagonale
typedef boost::numeric::ublas::triangular_matrix	<tipo_matrici, boost::numeric::ublas::upper >		MATRICI_TRIANGOLARI_UP;
typedef std::vector									<MATRICI_TRIANGOLARI_UP>							V_MATRICI_TRIANGOLARI_UP;
typedef std::vector									<V_MATRICI_TRIANGOLARI_UP>							V_V_MATRICI_TRIANGOLARI_UP;


//generali --> dense --> USARE CON ATTENZIONE
typedef boost::numeric::ublas::matrix				<tipo_matrici>										MATRICI;
typedef std::vector									<MATRICI>											V_MATRICI;
typedef std::vector									<V_MATRICI>											V_V_MATRICI;
typedef boost::numeric::ublas::matrix				<int>												MATRICI_INTERE;


//    mapped_matrix<double> m (3, 3, 3 * 3)
//							(righe, colonne, numero di elementi non nulli)
//boost::numeric::ublas::unbounded_array<double, std::allocator<double>
//typedef boost::numeric::ublas::compressed_matrix		<tipo_matrici>										MATRICI_SPARSE;

typedef boost::numeric::ublas::compressed_matrix		<tipo_matrici,boost::numeric::ublas::row_major,false,boost::numeric::ublas::unbounded_array<int> >  MATRICI_SPARSE;

typedef std::vector										<MATRICI_SPARSE>                                    V_MATRICI_SPARSE;
typedef std::vector										<V_MATRICI_SPARSE>									V_V_MATRICI_SPARSE;
typedef std::vector										<V_V_MATRICI_SPARSE>								V_V_V_MATRICI_SPARSE;

//si chiama con MATRICI_DIAGONALI(righe,colonne)
// (size_type size1, size_type size2, size_type lower = 0, size_type upper = 0)
typedef boost::numeric::ublas::banded_matrix		<tipo_matrici>										MATRICI_DIAGONALI;
typedef std::vector									<MATRICI_DIAGONALI>									V_MATRICI_DIAGONALI;
typedef std::vector									<V_MATRICI_DIAGONALI>								V_V_MATRICI_DIAGONALI;

typedef boost::numeric::ublas::identity_matrix		<tipo_matrici>										MATRICI_IDENTICHE;
typedef std::vector									<MATRICI_IDENTICHE>									V_MATRICI_IDENTICHE;
typedef std::vector									<V_MATRICI_IDENTICHE>								V_V_MATRICI_IDENTICHE;

//typedef boost::numeric::ublas::matrix				<tipo_matrici,boost::numeric::ublas::column_major>	MATRICI_CT;

typedef boost::numeric::ublas::vector				<tipo_matrici>						VETTORI;
//typedef boost::numeric::ublas::sparse_vector		<tipo_matrici>						VETTORI;
typedef std::vector									<VETTORI>							V_VETTORI;
typedef std::vector									<V_VETTORI>							V_V_VETTORI;

namespace ublas = boost::numeric::ublas;

#endif
