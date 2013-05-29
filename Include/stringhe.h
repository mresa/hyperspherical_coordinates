/*
 *  stringhe.h
 *
 *  Converte un numero in stringa, per poterlo usare nel titolo di un file
 *  Created by marco on 24/05/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef STRINGHE_H
#define STRINGHE_H

#include <sstream>
#include <string>
#include <stdexcept>

//classe di controllo per la conversione

class CattivaConversione : public std::runtime_error {
public:
	CattivaConversione(const std::string& s)
	: std::runtime_error(s)
	{ }
};
 
//funzione che scrive un numero in una stringa

inline std::string scrivi(const int & x) {
	std::ostringstream o;
	if (!(o << x))
		throw CattivaConversione("scrivi(int)");
	return o.str();
 }
 
// si usa come una stringa:
 
/* int main (int argc, char * const argv[]) {
	int n = 100;
    std::string s = "il valore è " + scrivi(n);
	
	std::ofstream myfile;

	myfile.open (s.c_str());
	myfile.close();
	return 0;
}
*/

#endif
