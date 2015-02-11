#include <RcppArmadillo.h>

using namespace arma;

cube prodint(cube & dna, int nstate, int ltimes) {

    cube aj(dna.begin(), nstate, nstate, ltimes);
    mat I = eye<mat>(nstate, nstate);

    aj.slice(0) = aj.slice(0) + I;

    for (int i = 1; i < ltimes; ++i) {
	aj.slice(i) = aj.slice(i-1) * (I + aj.slice(i));
    }

    return aj;
}

    
    

    
