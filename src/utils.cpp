#include <RcppArmadillo.h>

using namespace arma;

cube prodint(const cube & dna, int nstate, int ltimes) {

    cube aj(dna.begin(), nstate, nstate, ltimes);
    mat I = eye<mat>(nstate, nstate);

    aj.slice(0) = aj.slice(0) + I;

    for (int i = 1; i < ltimes; ++i) {
	aj.slice(i) = aj.slice(i-1) * (I + aj.slice(i));
    }

    return aj;
}

// cube prodint(const cube & dna, int nstate, int ltimes) {

//     cube aj(dna.begin(), nstate, nstate, ltimes);
//     mat I = eye<mat>(nstate, nstate);

//     aj.slice(0) = aj.slice(0) + I;

//     cube::slice_iterator c = aj.begin_slice(1);
//     cube::slice_iterator d = aj.end();
//     cube::iterator a = aj.begin();
//     cube::slice_iterator b = aj.begin_slice(ltimes - 1);

//     cube::iterator i; cube::iterator j;
//     for (i=c, j=a; i!=d; i += nstate * nstate, j += nstate * nstate) {
	
// 	mat bouh(i, nstate, nstate, false);
// 	mat yah(j, nstate, nstate, false);

// 	bouh = yah * (I + bouh);
//     }

//     return aj;
// }

cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes) {

    cube dna(nstate, nstate, ltimes);
    dna.zeros();
    
    for (int t=0; t < ltimes; ++t) {
	for (int i=0; i < nstate; ++i) {
	    if (nrisk.at(t, i) != 0) {
		for (int j=0; j< nstate; ++j) {
		    dna.at(i, j, t) = nev.at(i, j, t) / nrisk.at(t, i);
		}
	    }
	}
	
	mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	vec d = sum(tmp, 1);
    	tmp.diag() = -d;

    }
    
    return dna;

}
