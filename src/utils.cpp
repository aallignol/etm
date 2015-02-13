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

cube deltaNA(const cube & nev, const mat & nrisk, int nstate, int ltimes) {

    cube dna(nstate, nstate, ltimes);
    dna.zeros();
    
    for (int t=0; t < ltimes; ++t) {
	for (int i=0; i < nstate; ++i) {
	    for (int j=0; j< nstate; ++j) {
		if (nrisk.at(i, t) != 0) {
		    dna.at(i, j, t) = nev.at(i, j, t) / nrisk.at(i, t);
		}
	    }
	}
	
	mat tmp(dna.slice(t).begin(), nstate, nstate, false);
	vec d = sum(tmp, 1);
    	tmp.diag() = -d;

    }
    
    return dna;

}
