/*
 *  Test_rw_matrici.cpp
 *  Tesi
 *
 *  Created by Marco Resa on 13/05/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string.h>
#include <vector>
#include <fstream>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include "tipi_definiti.h"
#include "stringhe.h"
#include "io_sparse.h"


int main(int argc, char *argv[])
{
	using namespace boost::numeric::ublas;
	using std::cin;
	using std::cout;
	using std::endl;
	std::fstream f;
	std::string Nome_file_txt = "testfile.txt";
	std::string Nome_file_bin = "testfile.bin";
	const size_t size = 5;
	const size_t i_index[5] = { 0, 1, 2, 3, 4 };
	const size_t j_index[5] = { 0, 2, 0, 4, 4 };
	
	{
		MATRICI_SPARSE SM(size,size);
		for (size_t i=0; i<size; ++i) SM.insert_element(i_index[i], j_index[i], 1.0);
		f.open(Nome_file_txt.c_str(),std::ios::out);
		f << io::sparse(SM) << std::endl;
		f.close();
		cout << "Scritta matrice sparsa = " << SM << endl;
	}
	{
		MATRICI_SPARSE SM;
		f.open(Nome_file_txt.c_str(),std::ios::in);
		f >> io::sparse(SM);
		f.close();
		cout << "Letta matrice sparsa = " << SM << endl;
	}
	{
		MATRICI_SPARSE SM(size,size);
		for (size_t i=0; i<size; ++i) SM.insert_element(i_index[i], j_index[i], 1.0);
		f.open(Nome_file_bin.c_str(),std::ios::out | std::ios::binary);
		f.write(reinterpret_cast<char *>(&SM),sizeof(SM));
		f.close();
		cout << "Scritta binaria matrice sparsa = " << SM << endl;
	}
	{
		MATRICI_SPARSE SM;
		f.open(Nome_file_bin.c_str(),std::ios::in | std::ios::binary);
		f.read(reinterpret_cast<char *>(&SM),sizeof(SM));
		f.close();
		cout << "Letta binaria matrice sparsa = " << SM << endl;
	}
	
	return EXIT_SUCCESS;
}