/*
 *  Test_rw_vettori.cpp
 *  Tesi
 *
 *  Created by Marco Resa on 22/05/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

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
#include <boost/random.hpp>


int main(int argc, char *argv[])
{
	using namespace boost::numeric::ublas;
	using std::cin;
	using std::cout;
	using std::endl;
	std::fstream f;
	std::string Nome_file_txt = "testfile.txt";
	std::string Nome_file_bin = "testfile.bin";
	const size_t size = 100;
	typedef boost::lagged_fibonacci607 Gen_i;
	Gen_i mygen_i;
	
	{
		VETTORI Vec(size);
		std::generate(Vec.begin(),Vec.end(),mygen_i);
		f.open(Nome_file_txt.c_str(),std::ios::out);
		f << Vec << std::endl;
		f.close();
		cout << "Scritto vettore = " << Vec << endl;
	}
	{
		VETTORI Vec;
		cout << "NON Letto vettore = " << Vec << endl;
		f.open(Nome_file_txt.c_str(),std::ios::in);
		f >> Vec;
		f.close();
		cout << "Letto vettore = " << Vec << endl;
	}
/*	{
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
*/	
	return EXIT_SUCCESS;
}
